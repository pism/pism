import numpy as np

import PISM
import PISM.testing

ctx = PISM.Context()

def create_mask(shape):
    """Create the input mask and corresponding desired results."""
    assert shape[0] >= 11 and shape[1] >= 11

    mask = np.zeros(shape, dtype=int)
    mask[:] = PISM.MASK_ICE_FREE_OCEAN
    mask[6:,:4] = PISM.MASK_GROUNDED
    # not an iceberg
    mask[2:6,:2] = PISM.MASK_FLOATING
    mask_no_icebergs_fem = mask.copy()
    # add a patch that is not an iceberg int the FD context and *is* an iceberg in the FEM
    # context
    mask[:6,3:7] = PISM.MASK_FLOATING
    mask_no_icebergs_fd = mask.copy()
    # add a couple of icebergs (in both meanings)
    mask[6:,7:] = PISM.MASK_FLOATING
    mask[:3, 8:] = PISM.MASK_FLOATING

    return mask, mask_no_icebergs_fd, mask_no_icebergs_fem

def check(version):
    """Check that an iceberg remover correctly identifies and removes 'icebergs'."""

    grid = PISM.testing.shallow_grid(11, 11, Lx=1e4, Ly=1e4)

    cell_type     = PISM.IceModelVec2CellType(grid, "cell_type", PISM.WITH_GHOSTS)
    vel_bc_mask       = PISM.IceModelVec2Int(grid, "vel_bc_mask", PISM.WITH_GHOSTS)
    ice_thickness = PISM.IceModelVec2S(grid, "thk", PISM.WITHOUT_GHOSTS)

    input_mask, result_fd, result_fem = create_mask(cell_type.shape())

    with PISM.vec.Access(cell_type, ice_thickness):
        for i,j in grid.points():
            cell_type[i, j] = float(input_mask[j, i])
            ice_thickness[i, j] = cell_type.icy(i, j)

    cell_type.update_ghosts()

    if version == "fd":
        model = PISM.IcebergRemover(grid)
        desired_result = result_fd
    elif version == "fem":
        model = PISM.IcebergRemoverFEM(grid)
        desired_result = result_fem
    else:
        raise ValueError("invalid iceberg remover type: {}".format(version))

    model.update(vel_bc_mask, cell_type, ice_thickness)
    np.testing.assert_equal(cell_type.numpy(), desired_result)

def iceberg_remover_fd_test():
    """Iceberg remover (FD version)"""

    check("fd")

def iceberg_remover_fem_test():
    """Iceberg remover (FEM version)"""

    check("fem")
