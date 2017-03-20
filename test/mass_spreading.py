import PISM

log = PISM.Context().log

def disc(thickness, x0, y0, H, R):
    """Set ice thickness to H within the disc centered at (x0,y0) of
    radius R and 0 elsewhere.
    """

    grid = thickness.get_grid()

    with PISM.vec.Access(nocomm=thickness):
        for (i, j) in grid.points():
            x  = grid.x(i)
            y  = grid.y(j)
            d2 = (x - x0)**2 + (y - y0)**2
            r2 = R**2
            if d2 <= r2:
                thickness[i, j] = H
            else:
                thickness[i, j] = 0.0

    thickness.update_ghosts()

def set_ice_thickness(output):
    """Exact solution at time time. Corresponds to a disc that is rotated
    around the origin (one revolution per time unit)."""

    grid = output.get_grid()

    L = min(grid.Lx(), grid.Ly())

    R = 0.25 * L

    disc(output, 0, 0, 1, R)

def set_velocity(v):
    """Initialize the velocity field to a rigid rotation around the
    origin. This is slow, but it works.

    """
    grid = v.get_grid()

    scalar_velocity = 0.7

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

    Mx = 101
    My = 101

    grid = PISM.IceGrid_Shallow(ctx, 1, 1, 0, 0, Mx, My, PISM.NOT_PERIODIC)

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

    # grid info
    geometry.cell_area().set(grid.dx() * grid.dy())
    geometry.latitude().set(0.0)
    geometry.longitude().set(0.0)
    # environment
    geometry.bed_elevation().set(-10.0)
    geometry.sea_level_elevation().set(0.0)
    # ice
    set_ice_thickness(geometry.ice_thickness())
    geometry.ice_area_specific_volume().set(0.0)

    geometry.ensure_consistency(0.0)

    set_velocity(v)
    v_bc_mask.set(0.0)          # all points are B.C. points, but it does not matter
    set_ice_thickness(H_bc_mask)
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
        geometry.ice_thickness().dump("thk-%05d.nc" % (j+1))
        geometry.ice_area_specific_volume().dump("Href-%05d.nc" % (j+1))
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
        profiling.end("modify")

        t += dt
        j += 1
    profiling.stage_end("ge");

    profiling.report("profiling.py")

    return geometry, dt
