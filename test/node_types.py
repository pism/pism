import PISM
import time
import numpy as np
import pylab as plt

# This script tests compute_node_types() using a small 11*11 grid.
# Note that icy "tongues" of width 1 cannot be resolved by the Q1 FEM
# grid and so nodes corresponding to these are considered "exterior".

np.set_printoptions(precision=5, suppress=True, linewidth=100)

ctx = PISM.Context()

ice_density = ctx.config.get_double("constants.ice.density")
ocean_density = ctx.config.get_double("constants.sea_water.density")

mu = ice_density / ocean_density


def allocate_grid(ctx):
    params = PISM.GridParameters(ctx.config)
    params.Lx = 1e5
    params.Ly = 1e5
    params.Lz = 1000
    params.Mx = 11
    params.My = 11
    params.Mz = 5
    params.registration = PISM.CELL_CORNER
    params.periodicity = PISM.NOT_PERIODIC
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)


def allocate_storage(grid):
    ice_thickness = PISM.model.createIceThicknessVec(grid)

    mask = PISM.model.createIceMaskVec(grid)
    mask.set_name("node_type")

    return ice_thickness, mask


def print_vec(vec):
    v0 = vec.allocate_proc0_copy()
    vec.put_on_proc0(v0.get())

    shape = vec.get_dm().get().sizes

    print(vec.get_name())
    print(v0.get()[:].reshape(shape, order="f"))


def spy_vec(vec, value):
    v0 = vec.allocate_proc0_copy()
    vec.put_on_proc0(v0.get())

    shape = vec.get_dm().get().sizes

    array = v0.get()[:].reshape(shape, order="f")

    print(vec.get_name())
    print(array)

    plt.title(vec.get_name())
    plt.imshow(array, interpolation="nearest")


def init(H, vec):

    grid = vec.get_grid()

    K = 5
    R = 2

    with PISM.vec.Access(nocomm=[vec]):
        for (i, j) in grid.points():
            if abs(i - K) < R or abs(j - K) < R:
                vec[i, j] = H
            else:
                vec[i, j] = 0.0

            if abs(i - K) > R or abs(j - K) > R:
                vec[i, j] = 0.0

            if abs(i - K) < 1 or abs(j - K) < 1:
                vec[i, j] = H

            if abs(i - K) > R + 2 or abs(j - K) > R + 2:
                vec[i, j] = 0.0

    vec.update_ghosts()


def node_type_test():

    # allocation
    grid = allocate_grid(ctx)

    ice_thickness, mask = allocate_storage(grid)

    H = 1.0
    thickness_threshold = 0.5

    # initialization
    sea_level = 500.0

    init(H, ice_thickness)

    PISM.compute_node_types(ice_thickness, thickness_threshold, mask)

    # inspection / comparison
    plt.figure()
    spy_vec(ice_thickness, H)
    plt.figure()
    spy_vec(mask, 1.0)
    plt.show()


node_type_test()
