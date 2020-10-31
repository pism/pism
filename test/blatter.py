#!/usr/bin/env python3
"""This script runs the Blatter stress balance solver.
"""

import numpy as np
import PISM
import PISM.testing

ctx = PISM.Context()
config = ctx.config

config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")

def H(H_max, R_max, r):
    SperA = 31556926.0
    n = 3.0
    H0 = 3600.0
    R0 = 750000.0
    t = 450 * SperA

    # alpha=(2-(n+1)*lambda)/(5*n+3)=1/9
    alpha = 1.0 / 9.0
    # beta=(1+(2*n+1)*lambda)/(5*n+3)=1/18
    beta = 1.0 / 18.0
    # t0 = (beta/Gamma) * pow((2n+1)/(n+1),n)*(pow(R0,n+1)/pow(H0,2n+1))
    t0 = 422.45 * SperA

    Rmargin = R0 * pow(t / t0, beta);
    if (r < Rmargin):
        H = H0 * pow(t / t0, -alpha) * pow(1.0 - pow(pow(t / t0, -beta) * (r / R0), (n + 1) / n),
                                           n / (2*n + 1))
        return H
    else:
        return 0.0

def H0(H_max, R_max, r):
    return H_max * np.sqrt(max(1.0 - (r / R_max)**2, 0.0))

def allocate(Mx, Mz):
    H_max = 1000.0
    R_max = 750e3

    P = PISM.GridParameters(config)

    P.Mx = Mx
    P.Lx = 800e3
    P.x0 = 0.0

    dx = (2 * P.Lx) / (P.Mx - 1)

    P.Ly = dx
    P.My = 3
    P.y0 = 0.0

    P.registration = PISM.CELL_CORNER
    P.periodicity = PISM.Y_PERIODIC
    P.z = PISM.DoubleVector([0.0, 5000.0])

    P.ownership_ranges_from_options(ctx.size)

    grid = PISM.IceGrid(ctx.ctx, P)

    geometry = PISM.Geometry(grid)

    enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
    enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

    yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
    yield_stress.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

    with PISM.vec.Access(nocomm=[geometry.bed_elevation, geometry.ice_thickness]):
        for (i, j) in grid.points():
            r = abs(grid.x(i))
            geometry.bed_elevation[i, j] = 0.0
            geometry.ice_thickness[i, j] = H(H_max, R_max, r)

    geometry.sea_level_elevation.set(0.0)

    geometry.ensure_consistency(0.0)

    # this value is not important (we use an isothermal flow law)
    enthalpy.set(1e5)

    # this has to be high enough to prevent sliding
    yield_stress.set(10 * config.get_number("basal_yield_stress.constant.value"))

    return grid, geometry, enthalpy, yield_stress

def run(Mx, Mz, glen_exponent):
    grid, geometry, enthalpy, yield_stress = allocate(Mx, Mz)

    n_levels          = int(config.get_number("stress_balance.blatter.n_levels"))
    coarsening_factor = int(config.get_number("stress_balance.blatter.coarsening_factor"))

    config.set_number("stress_balance.blatter.Glen_exponent", glen_exponent)

    model = PISM.Blatter(grid, Mz, n_levels, coarsening_factor)

    inputs = PISM.StressBalanceInputs()

    inputs.geometry = geometry
    inputs.basal_yield_stress = yield_stress
    inputs.enthalpy = enthalpy

    model.update(inputs, True)

    return grid.x(), geometry.ice_thickness.numpy(), model.velocity_u_sigma().numpy()

if __name__ == "__main__":

    Mx = int(config.get_number("grid.Mx"))
    Mz = int(config.get_number("stress_balance.blatter.Mz"))

    x, u = run(Mx, Mz, 3)

    print(u.shape)

# check if shapes of elements play a role in this by changing the bed topography to match
# surface elevation (i.e. use constant ice thickness)

# try changing ice hardness and see what effect it has
