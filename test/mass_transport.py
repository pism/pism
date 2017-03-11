import numpy as np
import PISM

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
    x0, y0 = 0.5 * L, 0.0

    R = 0.25 * L
    x,y = M * np.matrix([x0, y0]).T

    disc(output, x, y, 1, R)

def set_velocity(v):
    """Initialize the velocity field to a rigid rotation around the
    origin. This is slow, but it works.

    """
    grid = v.get_grid()

    with PISM.vec.Access(nocomm=v):
        for (i, j) in grid.points():
            x  = grid.x(i)
            y  = grid.y(j)

            R = PISM.radius(grid, i, j)

            phi = np.arctan2(y, x)

            v[i, j].u = -R * np.sin(phi) * 2 * np.pi
            v[i, j].v =  R * np.cos(phi) * 2 * np.pi

def mass_transport_test(t_final, C=1.0):
    "Test GeometryEvolution::step()"

    ctx = PISM.Context().ctx

    Mx = 201
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
    geometry.bed_elevation().set(0.0)
    geometry.sea_level_elevation().set(0.0)
    # ice
    set_ice_thickness(geometry.ice_thickness(), time=0)
    geometry.ice_area_specific_volume().set(0.0)

    geometry.ensure_consistency(0.0)

    set_velocity(v)
    v_bc_mask.set(0.0)          # all points are B.C. points, but it does not matter
    H_bc_mask.set(0.0)
    SMB.set(0.0)
    BMR.set(0.0)

    geometry.ice_thickness().dump("thk-0.nc")

    t = 0.0
    dt = PISM.max_timestep_cfl_2d(geometry.ice_thickness(),
                                  geometry.cell_type(),
                                  v).dt_max.value() * C

    j = 0
    while t < t_final:
        print t

        ge.step(geometry, dt,
                v,
                Q,
                v_bc_mask,
                H_bc_mask,
                SMB,
                BMR)

        geometry.ice_thickness().add(1.0, ge.thickness_change_due_to_flow())

        geometry.ice_thickness().dump("thk-%05d.nc" % j)

        if t + dt > t_final:
            dt = t_final - t

        t += dt
        j += 1

    geometry.ice_thickness().dump("thk-1.nc")

    return ge, dt

def move_particle(velocity, t_final, C=1.0):
    """Update position of a particle using the provided velocity field.
    Used to make sure we really do have a rotational velocity field.
    """
    grid = velocity.get_grid()
    L = min(grid.Lx(), grid.Ly())

    # compute max time step
    H = PISM.model.createIceThicknessVec(grid)
    m = PISM.model.createIceMaskVec(grid)

    H.set(100.0)
    m.set(PISM.MASK_GROUNDED)

    dt = PISM.max_timestep_cfl_2d(H, m, velocity).dt_max.value() * C

    # initial position
    p0 = PISM.Vector2(L / 2.0, 0.0)

    p = PISM.Vector2(p0.u, p0.v)
    t = 0.0
    with PISM.vec.Access(nocomm=velocity):
        while t < t_final:
            if t + dt > t_final:
                dt = t_final - t

            v = velocity.interpolate(p.u, p.v)

            p += v * dt

            t += dt

    return (p - p0).magnitude()
