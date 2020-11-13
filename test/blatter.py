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

def H1(H_max, R_max, r):
    n = 4.0 / 3.0
    return H_max * np.sqrt(max(1.0 - (r / R_max)**n, 0.0))

def H2(H_max, R_max, r):
    if r < 0.5 * R_max:
        return 0.5 * H_max
    else:
        return H_max * max(1.0 - r / R_max, 0.0)

def allocate(Mx, Mz):
    H_max = 500.0
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

# try
#
# - refinement in x
#
# - refinement in z
#
# - changing ice hardness
#
# - changing the Glen exponent
#
# - changing the viscosity regularization parameter
#
# - changing the basal yield stress
#
# - extending the ice at grounded margins by an extra (very thin) "buffer", moving the
#   margin by 1 grid point and reducing the thickness of the ice at the margin (sort of
#   like the implementation in MALI)
#
# - changing bed topography and ice thickness (keeping surface elevation constant) so that
#   all elements have approximately the same shape (right now elements at the margin are
#   more deformed than the ones near the dome)
#
# - solving the same problem, but with the surface gradient computed *analytically* at
#   quadrature points
#
# - changing the number of quadrature points
#
# - Observing the maximum dh/dx as the grid is refined in x. Is there a relationship
#   between max. slope and the presence (absence) of the spike?
#
# - it looks like element distortion plays a role. I need to check if ice thicknesses such
#   that H_{n+1} <= 0.75 H_n *in the middle of the domain* cause issues. This would
#   indicate that lateral BC (and the accuracy of prescribing it) does not explain the
#   spike in ice speed at the margin.
#
# - I need to make sure that the lateral BC is correct: the BP model has the cryostatic
#   pressure built in, so I'm not sure we need to have it in the BC. Gotta check this.
#
# - I need to add a flag that toggles using the lateral BC at grounded margins.
#
# - I need to set up an X-Z test using the van der Veen profile. In this case I have exact
#   SSA velocities corresponding to the (known to be correct) lateral BC for the SSA
#   system. The BP solution should be close to the SSA. This should tell me if the BP
#   lateral BC I have are wrong (in the sense of "using the right equations", not
#   "correctly implementing chosen equations").
