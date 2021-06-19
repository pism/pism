#!/usr/bin/env python3
"""This script runs ISMIP-HOM experiments A-D using the Blatter stress balance solver.

For experiments A and C and length scales of 5 and 10 km it is helpful to use a good
preconditioner. Multigrid using semi-coarsening in the vertical direction appears to be
very effective. At higher grid resolutions (i.e. when the coarse problem still has a large
number of unknowns) using algebraic multigrid preconditioner on the coarse level may help.

Here's an example:

mpiexec -n 2 python3 run-ismiphom.py \
   -pc_type mg \
   -pc_mg_levels 2 \
   -stress_balance.blatter.coarsening_factor 4 \
   -snes_monitor \
   -ksp_monitor

This command uses "aggressive" coarsening (factor of 4) and 2 multigrid levels (5 and 2
nodes in the vertical direction).

To use AMG on the coarse level, add "-mg_coarse_pc_type gamg".

"""

import numpy as np
import PISM
import PISM.testing
from PISM.util import convert

seconds_per_year = convert(1.0, "year", "second")

ctx = PISM.Context()
config = ctx.config

tests = {"A" : PISM.HOM_A,
         "B" : PISM.HOM_B,
         "C" : PISM.HOM_C,
         "D" : PISM.HOM_D}

def surface_AB(x, y, L):
    alpha = 0.5 * (np.pi / 180.0) # 0.5 degrees
    return -x * np.tan(alpha)

def surface_C(x, y, L):
    alpha = 0.1 * (np.pi / 180.0) # 0.1 degrees
    return -x * np.tan(alpha)

def bed_A(x, y, L):
    omega = 2 * np.pi / L
    return surface_AB(x, y, L) - 1000.0 + 500.0 * np.sin(omega * x) * np.sin(omega * y)

def bed_B(x, y, L):
    omega = 2 * np.pi / L
    return surface_AB(x, y, L) - 1000.0 + 500.0 * np.sin(omega * x)

def bed_C(x, y, L):
    return surface_C(x, y, L) - 1000.0

def tauc_AB(x, y, L):
    return 1e16             # approximates infinity

def tauc_C(x, y, L):
    omega = 2.0 * np.pi / L
    return 1000.0 + 1000 * np.sin(omega * x) * np.sin(omega * y)

def tauc_D(x, y, L):
    omega = 2.0 * np.pi / L
    return 1000.0 + 1000 * np.sin(omega * x)

s =  {"A" : surface_AB,
      "B" : surface_AB,
      "C" : surface_C,
      "D" : surface_C}

b =  {"A" : bed_A,
      "B" : bed_B,
      "C" : bed_C,
      "D" : bed_C}

tauc = {"A" : tauc_AB,
        "B" : tauc_AB,
        "C" : tauc_C,
        "D" : tauc_D}

def init(testname, L):

    # the domain is a square [0, L]*[0, L]
    P = PISM.GridParameters(config)
    P.Lx = L / 2.0
    P.Ly = L / 2.0
    P.x0 = L / 2.0
    P.y0 = L / 2.0
    P.periodicity = PISM.XY_PERIODIC
    P.z = PISM.DoubleVector(np.linspace(0, 2000, 201))
    P.horizontal_size_from_options()

    # use 3 grid points in the y direction for x-z tests
    if testname in "BD":
        P.My = 3

    P.ownership_ranges_from_options(ctx.size)

    grid = PISM.IceGrid(ctx.ctx, P)

    geometry = PISM.Geometry(grid)

    grid.variables().add(geometry.ice_thickness)
    grid.variables().add(geometry.cell_type)
    grid.variables().add(geometry.bed_elevation)

    # enthalpy values are irrelevant: these tests use an isothermal flow law
    enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
    enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

    yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
    yield_stress.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

    with PISM.vec.Access([yield_stress, geometry.ice_thickness, geometry.bed_elevation]):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            # Convert from "Pa year / m" to "Pa s / m"
            yield_stress[i, j] = tauc[testname](x, y, L) * seconds_per_year
            # Set ice thickness: it is not used by the Blatter-Pattyn solver itself but is
            # needed to compute vertically-averaged ice velocity, etc.
            geometry.ice_thickness[i, j] = s[testname](x, y, L) - b[testname](x, y, L)
            geometry.bed_elevation[i, j] = b[testname](x, y, L)

    geometry.sea_level_elevation.copy_from(geometry.bed_elevation)
    # ensure that the whole domain is grounded
    geometry.sea_level_elevation.shift(-100.0)

    geometry.ensure_consistency(0.0)

    return grid, geometry, enthalpy, yield_stress

def set_constants(config):

    # Set up the sliding law so that tauc == beta**2:
    config.set_flag("basal_resistance.pseudo_plastic.enabled", True)
    config.set_number("basal_resistance.pseudo_plastic.q", 1.0)
    config.set_number("basal_resistance.pseudo_plastic.u_threshold",
                      convert(1.0, "m / s", "m / year"))

    # Set flow law parameters
    config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")
    config.set_number("stress_balance.blatter.Glen_exponent", 3.0)
    config.set_number("flow_law.isothermal_Glen.ice_softness",
                      convert(1e-16, "Pa-3 year-1", "Pa-3 s-1"))

    # Set constants:
    config.set_number("constants.ice.density", 910.0)
    config.set_number("constants.standard_gravity", 9.81)

    # Support compression
    config.set_string("output.format", "netcdf4_serial")

def run_test(test_name, L, output_file):
    """Run an ISMIP-HOM test on a [0,L]*[0,L] domain and save results to an output_file."""

    grid, geometry, enthalpy, yield_stress = init(test_name, L)

    Mz = int(config.get_number("stress_balance.blatter.Mz"))

    coarsening_factor = int(config.get_number("stress_balance.blatter.coarsening_factor"))

    model = PISM.BlatterISMIPHOM(grid, Mz, coarsening_factor, tests[test_name])

    stress_balance = PISM.StressBalance(grid, model, PISM.BlatterMod(model))

    stress_balance.init()

    inputs = PISM.StressBalanceInputs()

    inputs.geometry = geometry
    inputs.basal_yield_stress = yield_stress
    inputs.enthalpy = enthalpy

    stress_balance.update(inputs, True)

    try:
        output = PISM.util.prepare_output(output_file)
        output.set_compression_level(1)

        ds = stress_balance.diagnostics()
        ds["velsurf"].compute().write(output)

        geometry.ice_thickness.write(output)
        geometry.bed_elevation.write(output)
        geometry.ice_surface_elevation.write(output)

        yield_stress.write(output)
    finally:
        output.close()

if __name__ == "__main__":
    set_constants(config)

    for test in ["A", "B", "C", "D"]:
        for L in [5, 10, 20, 40, 80, 160]:
            run_test(test, 1e3 * L, "pism-hom-{}-{:03d}.nc".format(test, int(L)))
