import numpy as np
import PISM

"""This is a half-baked test for the mass transport code that sets up
a rotational velocity field and then steps forward in time for one revolution.

It is not ready to be used as a regression or a verification test. All
I can see right now is that our mass transport scheme is very
diffusive, which is not a surprise."""

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

def set_ice_thickness(output, time):
    """Exact solution at time time. Corresponds to a disc that is rotated
    around the origin (one revolution per time unit)."""

    grid = output.get_grid()

    L = min(grid.Lx(), grid.Ly())

    # 1 revolution per 1 time unit
    phi = 2 * np.pi * time

    M = np.matrix([[np.cos(phi), -np.sin(phi)],
                   [np.sin(phi),  np.cos(phi)]])

    # center coordinates at time 0
    x0, y0 = 0.25 * L, 0.0

    R = 0.25 * L
    x,y = M * np.matrix([x0, y0]).T

    disc(output, x, y, 1, R)

def set_velocity(v):
    """Initialize the velocity field to a rigid rotation around the
    origin.

    """
    grid = v.get_grid()

    radial_velocity = 2 * np.pi

    with PISM.vec.Access(nocomm=v):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)

            v[i, j].u = -y * radial_velocity
            v[i, j].v =  x * radial_velocity

    v.update_ghosts()

def mass_transport_test(t_final, C=1.0):
    "Test GeometryEvolution::step()"

    ctx = PISM.Context().ctx

    config = PISM.Context().config

    # config.set_boolean("geometry.part_grid.enabled", True)

    Mx = 101
    My = 101

    grid = PISM.IceGrid_Shallow(ctx, 1, 1, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

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
    set_ice_thickness(geometry.ice_thickness(), time=1)
    geometry.ice_area_specific_volume().set(0.0)

    geometry.ensure_consistency(0.0)

    set_velocity(v)
    v_bc_mask.set(0.0)          # all points are velocity B.C. points (but it does not matter)
    H_bc_mask.set(0.0)
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
        geometry.ensure_consistency(0.0)
        profiling.end("modify")

        t += dt
        j += 1
    profiling.stage_end("ge");

    profiling.report("profiling.py")

    return geometry, dt
