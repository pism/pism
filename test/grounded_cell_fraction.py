"""Check if the code computing the fraction of grid cells that is grounded produces
symmetric outputs for symmetric inputs.
"""

import PISM
import numpy as np

ctx = PISM.Context()

def init(mu, L, sea_level, thickness, test_type="box"):
    """Set ice thickness for the grounded fraction symmetry test.

    Parameters
    ----------
    mu :
        the ratio of densities of ice and ocean water
    L :
        controls the ice thickness near the center of the domain
    sea_level :
        the sea level elevation
    thickness :
        the output argument
    test_type :
        controls the shape of the icy patch ("box" or "cross")

    In the old implementation based on the ad hoc extension of the 1D LI parameterization
    of grounding line position L corresponded to the grounded fraction in a cell
    containing the grounding line. It does not have a clear meaning in the current
    context, other than 0 corresponds to floatation thickness and 1 to double that.

    """
    k = {0.0: 8,
         0.25: 7,
         0.5: 6,
         0.75: 5,
         1.0: 4}

    H0 = (8.0 / k[L]) * (sea_level / mu)
    H1 = 0.5 * H0

    grid = thickness.grid()

    with PISM.vec.Access(thickness):
        for (i, j) in grid.points():
            if test_type == "box" and abs(i - 3) < 2 and abs(j - 3) < 2:
                thickness[i, j] = H0
            elif test_type == "cross" and abs(i - 3) < 2 and abs(j - 3) < 2 and (i == 3 or j == 3):
                thickness[i, j] = H0
            else:
                thickness[i, j] = H1

            if abs(i - 3) >= 3 or abs(j - 3) >= 3:
                thickness[i, j] = 0.0

    thickness.update_ghosts()

def run(L, test):
    """Calls PISM's code to compute grounded cell fraction.

    Parameters
    ----------

    L :
        controls the thickness in the interior of the icy patch
    test :
        chooses the shape of the patch ("box" or "cross")
    """
    Lx, Ly = 1e5, 1e5
    Mx, My = 7, 7
    x0, y0 = 0.0, 0.0
    grid = PISM.Grid.Shallow(ctx.ctx,
                                Lx, Ly, x0, y0, Mx, My,
                                PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    sea_level = 500.0

    geometry = PISM.Geometry(grid)

    ice_density = ctx.config.get_number("constants.ice.density")
    ocean_density = ctx.config.get_number("constants.sea_water.density")

    init(ice_density / ocean_density, L, sea_level, geometry.ice_thickness, test)

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(sea_level)

    # this call computes the grounded fraction of each cell
    geometry.ensure_consistency(0.0)

    return geometry

def check_symmetry(var):
    """Check symmetry of a NumPy array"""
    np.testing.assert_almost_equal(var, np.flipud(var))
    np.testing.assert_almost_equal(var, np.fliplr(var))
    np.testing.assert_almost_equal(var, np.flipud(np.fliplr(var)))

def grounded_cell_fraction_test():
    """Check grounded cell fraction symmetry for symmetric inputs"""
    for test in ["box", "cross"]:
        for L in [0.0, 0.25, 0.5, 0.75, 1.0]:
            geometry = run(L, test)

            check_symmetry(geometry.cell_grounded_fraction.to_numpy())

if __name__ == "__main__":
    grounded_cell_fraction_test()
