import PISM
import math

log = PISM.Context().log

def disc(thickness, x0, y0, H0, R_inner, R_outer):
    """Set ice thickness to H0 within the disc centered at (x0,y0) of
    radius R_inner, C/R in an annulus R_inner < r <= R_outer and 0
    elsewhere.

    """

    grid = thickness.get_grid()

    R_inner_2 = R_inner**2
    R_outer_2 = R_outer**2

    C = H0 * R_inner

    with PISM.vec.Access(nocomm=thickness):
        for (i, j) in grid.points():
            x  = grid.x(i)
            y  = grid.y(j)
            d2 = (x - x0)**2 + (y - y0)**2
            if d2 <= R_inner_2:
                thickness[i, j] = H0
            elif d2 <= R_outer_2:
                thickness[i, j] = C / math.sqrt(d2)
            else:
                thickness[i, j] = 0.0

    thickness.update_ghosts()

def set_velocity(scalar_velocity, v):
    """Initialize the velocity field to a rigid rotation around the
    origin. This is slow, but it works.

    """
    grid = v.get_grid()

    with PISM.vec.Access(nocomm=v):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)

            r = max(PISM.radius(grid, i, j), 0.001)

            v[i, j].u = scalar_velocity * x / r
            v[i, j].v = scalar_velocity * y / r

    v.update_ghosts()

def mass_transport_test(t_final, C=1.0):
    "Test GeometryEvolution::step()"

    ctx = PISM.Context().ctx

    config = PISM.Context().config

    config.set_boolean("geometry.part_grid.enabled", True)

    Mx = 401
    My = 401

    grid = PISM.IceGrid_Shallow(ctx, 1, 1, 0, 0, Mx, My, PISM.NOT_PERIODIC)

    L = min(grid.Lx(), grid.Ly())

    R_inner = 0.25 * L

    geometry = PISM.Geometry(grid)

    v = PISM.model.create2dVelocityVec(grid)

    Q = PISM.IceModelVec2Stag()
    Q.create(grid, "Q", PISM.WITHOUT_GHOSTS)

    v_bc_mask = PISM.IceModelVec2Int()
    v_bc_mask.create(grid, "v_bc_mask", PISM.WITHOUT_GHOSTS)

    H_bc_mask = PISM.IceModelVec2Int()
    H_bc_mask.create(grid, "H_bc_mask", PISM.WITHOUT_GHOSTS)

    SMB = PISM.IceModelVec2S()
    SMB.create(grid, "SMB", PISM.WITHOUT_GHOSTS)

    BMR = PISM.IceModelVec2S()
    BMR.create(grid, "BMR", PISM.WITHOUT_GHOSTS)

    ge = PISM.GeometryEvolution(grid)

    spreading_velocity = 0.7

    # grid info
    geometry.cell_area().set(grid.dx() * grid.dy())
    geometry.latitude().set(0.0)
    geometry.longitude().set(0.0)
    # environment
    geometry.bed_elevation().set(-10.0)
    geometry.sea_level_elevation().set(0.0)
    # ice
    # save exact ice thickness
    R_outer = R_inner + spreading_velocity * t_final
    disc(geometry.ice_thickness(), 0, 0, 1, R_inner, R_outer)
    geometry.ice_thickness().dump("thk-exact.nc")
    # set initial ice thickness
    disc(geometry.ice_thickness(), 0, 0, 1, R_inner, R_inner)
    geometry.ice_area_specific_volume().set(0.0)

    geometry.ensure_consistency(0.0)

    set_velocity(spreading_velocity, v)
    v_bc_mask.set(0.0)          # all points are B.C. points, but it does not matter
    disc(H_bc_mask, 0, 0, 1, R_inner, R_inner)
    SMB.set(0.0)
    BMR.set(0.0)

    profiling = ctx.profiling()
    profiling.start()

    t = 0.0
    j = 0
    profiling.stage_begin("ge");
    while t < t_final:
        dt = PISM.max_timestep_cfl_2d(geometry.ice_thickness(),
                                      geometry.cell_type(),
                                      v).dt_max.value() * C

        if t + dt > t_final:
            dt = t_final - t

        log.message(2, "{}, {}\n".format(t, dt))

        profiling.begin("dump");
        # geometry.ice_thickness().dump("thk-%05d.nc" % (j+1))
        # geometry.ice_area_specific_volume().dump("Href-%05d.nc" % (j+1))
        profiling.end("dump")

        profiling.begin("step")
        ge.step(geometry, dt,
                v,
                Q,
                v_bc_mask,
                H_bc_mask,
                SMB,
                BMR)
        profiling.end("step")

        profiling.begin("modify")
        geometry.ice_thickness().add(1.0, ge.thickness_change_due_to_flow())
        geometry.ice_area_specific_volume().add(1.0, ge.area_specific_volume_change_due_to_flow())
        geometry.ensure_consistency(0.0)
        profiling.end("modify")

        t += dt
        j += 1
    profiling.stage_end("ge");

    # combine stuff
    geometry.ice_thickness().add(1.0, geometry.ice_area_specific_volume())
    # dump stuff
    geometry.ice_thickness().dump("thk-final.nc")

    profiling.report("profiling.py")

    return geometry, dt
