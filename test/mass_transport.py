import PISM
import numpy as np

""" Test mass transport code using a radially-symmetric setup in which
a disc expands uniformly in all directions. Given mass conservation,
the thickness of this disc outside of the fixed circular area (the
thickness Dirichlet B.C. area) is inversely proportional to the
distance from the center. Here we use this 'exact solution' to test
the symmetry of the produced ice thickness and linear convergence
towards the exact solution."""

log = PISM.Context().log


def disc(thickness, x0, y0, H0, R_inner, R_outer):
    """Set ice thickness to H0 within the disc centered at (x0,y0) of
    radius R_inner, C/R in an annulus R_inner < r <= R_outer and 0
    elsewhere.

    """

    grid = thickness.grid()

    R_inner_2 = R_inner**2
    R_outer_2 = R_outer**2

    C = H0 * R_inner

    with PISM.vec.Access(nocomm=thickness):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            d2 = (x - x0)**2 + (y - y0)**2
            if d2 <= R_inner_2:
                thickness[i, j] = H0
            elif d2 <= R_outer_2:
                thickness[i, j] = C / np.sqrt(d2)
            else:
                thickness[i, j] = 0.0

    thickness.update_ghosts()


def set_velocity(scalar_velocity, v):
    """Initialize the velocity field to a rigid rotation around the
    origin. This is slow, but it works.

    """
    grid = v.grid()

    with PISM.vec.Access(nocomm=v):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)

            r = max(PISM.radius(grid, i, j), 0.001)

            v[i, j].u = scalar_velocity * x / r
            v[i, j].v = scalar_velocity * y / r

    v.update_ghosts()


def run(Mx, My, t_final, part_grid, C=1.0):
    "Test GeometryEvolution::step()"

    ctx = PISM.Context().ctx

    config = PISM.Context().config

    config.set_boolean("geometry.part_grid.enabled", part_grid)

    grid = PISM.IceGrid_Shallow(ctx, 1, 1, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    assert t_final <= 1.0

    L = min(grid.Lx(), grid.Ly())
    R_inner = 0.25 * L
    spreading_velocity = 0.7
    R_outer = R_inner + spreading_velocity * t_final

    geometry = PISM.Geometry(grid)

    v         = PISM.IceModelVec2V(grid, "velocity", PISM.WITHOUT_GHOSTS)
    Q         = PISM.IceModelVec2Stag(grid, "Q", PISM.WITHOUT_GHOSTS)
    v_bc_mask = PISM.IceModelVec2Int(grid, "v_bc_mask", PISM.WITHOUT_GHOSTS)
    H_bc_mask = PISM.IceModelVec2Int(grid, "H_bc_mask", PISM.WITHOUT_GHOSTS)

    ge = PISM.GeometryEvolution(grid)

    # grid info
    geometry.latitude.set(0.0)
    geometry.longitude.set(0.0)
    # environment
    geometry.bed_elevation.set(-10.0)
    geometry.sea_level_elevation.set(0.0)
    # set initial ice thickness
    disc(geometry.ice_thickness, 0, 0, 1, R_inner, R_inner)
    geometry.ice_area_specific_volume.set(0.0)

    geometry.ensure_consistency(0.0)

    set_velocity(spreading_velocity, v)
    v_bc_mask.set(0.0)
    disc(H_bc_mask, 0, 0, 1, R_inner, R_inner)

    profiling = ctx.profiling()
    profiling.start()

    t = 0.0
    j = 0
    profiling.stage_begin("ge")
    while t < t_final:
        dt = PISM.max_timestep_cfl_2d(geometry.ice_thickness,
                                      geometry.cell_type,
                                      v).dt_max.value() * C

        if t + dt > t_final:
            dt = t_final - t

        log.message(2, "{}, {}\n".format(t, dt))

        profiling.begin("step")
        ge.flow_step(geometry, dt,
                     v,
                     Q,
                     v_bc_mask,
                     H_bc_mask)
        profiling.end("step")

        profiling.begin("modify")
        ge.apply_flux_divergence(geometry)
        geometry.ensure_consistency(0.0)
        profiling.end("modify")

        t += dt
        j += 1
    profiling.stage_end("ge")

    profiling.report("profiling_%d_%d.py" % (Mx, My))

    return geometry


def average_error(N):
    t_final = 1.0
    C = 1.0

    log.disable()
    geometry = run(N, N, t_final, True, C)
    log.enable()
    # combine stuff stored as thickness and as area specific volume
    geometry.ice_thickness.add(1.0, geometry.ice_area_specific_volume)

    grid = geometry.ice_thickness.grid()

    diff = PISM.IceModelVec2S(grid, "difference", PISM.WITHOUT_GHOSTS)
    exact = PISM.IceModelVec2S(grid, "thk", PISM.WITHOUT_GHOSTS)

    L = min(grid.Lx(), grid.Ly())
    R_inner = 0.25 * L
    spreading_velocity = 0.7
    R_outer = R_inner + spreading_velocity * t_final

    disc(exact, 0, 0, 1, R_inner, R_outer)

    exact.add(-1.0, geometry.ice_thickness, diff)

    # return the average error
    return diff.norm(PISM.PETSc.NormType.N1) / (N*N)


def part_grid_convergence_test():
    "Test that the error does go down as O(1/N)"

    np.testing.assert_almost_equal([average_error(N) for N in [51, 101]],
                                   [0.0338388,  0.0158498])


def part_grid_symmetry_test():
    """The initial condition and the velocity fields are radially
    symmetric, so the result should be too."""

    N = 51

    log.disable()
    geometry = run(N, N, 1, True, 1.0)
    log.enable()

    # combine stuff stored as thickness and as area specific volume
    geometry.ice_thickness.add(1.0, geometry.ice_area_specific_volume)

    # convert ice thickness to a NumPy array on rank 0 -- that way we can use flipud() and
    # fliplr().
    H = geometry.ice_thickness.numpy()

    np.testing.assert_almost_equal(H, np.flipud(H))
    np.testing.assert_almost_equal(H, np.fliplr(H))
    np.testing.assert_almost_equal(H, np.flipud(np.fliplr(H)))
