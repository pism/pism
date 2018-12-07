import PISM
import time
import numpy as np

np.set_printoptions(precision=5, suppress=True)

ctx = PISM.Context()

ice_density = ctx.config.get_double("constants.ice.density")
ocean_density = ctx.config.get_double("constants.sea_water.density")

mu = ice_density / ocean_density


def allocate_grid(ctx):
    params = PISM.GridParameters(ctx.config)
    params.Lx = 1e5
    params.Ly = 1e5
    params.Lz = 1000
    params.Mx = 7
    params.My = 7
    params.Mz = 5
    params.periodicity = PISM.NOT_PERIODIC
    params.registration = PISM.CELL_CORNER
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)


def allocate_storage(grid):
    ice_thickness = PISM.model.createIceThicknessVec(grid)

    # not used, but needed by GeometryCalculator::compute()
    surface = PISM.model.createIceSurfaceVec(grid)

    bed_topography = PISM.model.createBedrockElevationVec(grid)

    mask = PISM.model.createIceMaskVec(grid)

    gl_mask = PISM.model.createGroundingLineMask(grid)
    gl_mask_x = PISM.model.createGroundingLineMask(grid)
    gl_mask_x.set_name("gl_mask_x")
    gl_mask_y = PISM.model.createGroundingLineMask(grid)
    gl_mask_y.set_name("gl_mask_y")

    sea_level = PISM.model.createIceThicknessVec(grid)
    sea_level.set_name("sea_level")

    return ice_thickness, bed_topography, surface, mask, gl_mask, gl_mask_x, gl_mask_y, sea_level


def compute_mask(sea_level, bed_topography, ice_thickness, mask, surface):
    gc = PISM.GeometryCalculator(ctx.config)
    gc.compute(sea_level, bed_topography, ice_thickness, mask, surface)


def print_vec(vec):
    v0 = vec.allocate_proc0_copy()
    vec.put_on_proc0(v0.get())

    shape = vec.get_dm().get().sizes

    print(vec.get_name())
    print(v0.get()[:].reshape(shape, order="f"))


def init(mu, L, sea_level, vec, type="box"):
    k = {0.0: 8,
         0.25: 7,
         0.5: 6,
         0.75: 5,
         1.0: 4}

    H0 = (8.0 / k[L]) * (sea_level / mu)
    H1 = 0.5 * H0

    grid = vec.grid()

    with PISM.vec.Access(nocomm=[vec]):
        for (i, j) in grid.points():
            if type == "box" and abs(i - 3) < 2 and abs(j - 3) < 2:
                vec[i, j] = H0
            elif type == "cross" and abs(i - 3) < 2 and abs(j - 3) < 2 and (i == 3 or j == 3):
                vec[i, j] = H0
            else:
                vec[i, j] = H1

            if abs(i - 3) >= 3 or abs(j - 3) >= 3:
                vec[i, j] = 0.0

    vec.update_ghosts()


def grounded_cell_fraction_test():

    # allocation
    grid = allocate_grid(ctx)

    ice_thickness, bed_topography, surface, mask, gl_mask, gl_mask_x, gl_mask_y, _ = allocate_storage(grid)

    bed_topography.set(0.0)

    # initialization
    sea_level = 500.0
    for L in [0.0, 0.25, 0.5, 0.75, 1.0]:
        init(mu, L, sea_level, ice_thickness, "box")

        compute_mask(sea_level, bed_topography, ice_thickness, mask, surface)

        # computation of gl_mask
        PISM.compute_grounded_cell_fraction(ice_density, ocean_density, sea_level,
                                            ice_thickness, bed_topography, mask, gl_mask,
                                            gl_mask_x, gl_mask_y)

        # inspection / comparison
        print("L = %f" % L)
        print_vec(mask)
        print_vec(gl_mask_x)
        print_vec(gl_mask_y)
        print_vec(gl_mask)


def new_grounded_cell_fraction_test():

    # allocation
    grid = allocate_grid(ctx)

    ice_thickness, bed_topography, _, _, gl_mask, _, _, sea_level = allocate_storage(grid)

    # initialization
    bed_topography.set(0.0)
    sl = 500.0
    sea_level.set(sl)
    for L in [0.0, 0.25, 0.5, 0.75, 1.0]:
        init(mu, L, sl, ice_thickness, "box")

        # computation of gl_mask
        PISM.compute_grounded_cell_fraction(ice_density, ocean_density,
                                            sea_level,
                                            ice_thickness,
                                            bed_topography,
                                            gl_mask)

        # inspection / comparison
        print("L = %f" % L)
        print_vec(gl_mask)


grounded_cell_fraction_test()
new_grounded_cell_fraction_test()
