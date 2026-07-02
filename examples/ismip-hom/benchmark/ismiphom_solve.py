#!/usr/bin/env python3
"""Single ISMIP-HOM Blatter solve, instrumented for solver benchmarking.

Runs ONE ISMIP-HOM experiment as a single Blatter stress-balance solve and
reports the solve wall-clock time. Unlike a frozen-bed dome (e.g. eisII A),
ISMIP-HOM produces genuine 3-D velocities -- experiment C has spatially varying
basal sliding -- so the Blatter/BP linear solve (and its coarse-grid solve) is
actually exercised, which is the point of the benchmark.

All PETSc/PISM solver options (e.g. -bp_pc_type mg, -bp_mg_coarse_pc_type gamg)
are read from the command line; the grid size comes from -Mx/-My.

Controlled by environment variables (set by the benchmark driver):
  BENCH_TEST   ISMIP-HOM experiment: A, B, C, or D   (default C)
  BENCH_L_KM   domain side length, in km             (default 80)
  BENCH_MX     grid points in x                       (default 201)
  BENCH_MY     grid points in y                       (default 201)
  BENCH_OUT    output NetCDF file                     (default ismiphom_out.nc)

Adapted from examples/ismip-hom/abcd/run-ismiphom.py.
"""
import os
import time

import numpy as np
import PISM
from PISM.util import convert
from petsc4py import PETSc

seconds_per_year = convert(1.0, "year", "second")

TEST = os.environ.get("BENCH_TEST", "C").upper()
L = float(os.environ.get("BENCH_L_KM", "80")) * 1.0e3
MX = int(os.environ.get("BENCH_MX", "201"))
MY = int(os.environ.get("BENCH_MY", "201"))
OUT = os.environ.get("BENCH_OUT", "ismiphom_out.nc")

ctx = PISM.Context()
config = ctx.config

tests = {"A": PISM.HOM_A, "B": PISM.HOM_B, "C": PISM.HOM_C, "D": PISM.HOM_D}


def surface_AB(x, y, L):
    return -x * np.tan(0.5 * np.pi / 180.0)   # 0.5 degrees


def surface_C(x, y, L):
    return -x * np.tan(0.1 * np.pi / 180.0)   # 0.1 degrees


def bed_A(x, y, L):
    omega = 2 * np.pi / L
    return surface_AB(x, y, L) - 1000.0 + 500.0 * np.sin(omega * x) * np.sin(omega * y)


def bed_B(x, y, L):
    omega = 2 * np.pi / L
    return surface_AB(x, y, L) - 1000.0 + 500.0 * np.sin(omega * x)


def bed_C(x, y, L):
    return surface_C(x, y, L) - 1000.0


def tauc_AB(x, y, L):
    return 1e16             # approximates infinity (no sliding)


def tauc_C(x, y, L):
    omega = 2.0 * np.pi / L
    return 1000.0 + 1000.0 * np.sin(omega * x) * np.sin(omega * y)


def tauc_D(x, y, L):
    omega = 2.0 * np.pi / L
    return 1000.0 + 1000.0 * np.sin(omega * x)


surface = {"A": surface_AB, "B": surface_AB, "C": surface_C, "D": surface_C}
bed = {"A": bed_A, "B": bed_B, "C": bed_C, "D": bed_C}
tauc = {"A": tauc_AB, "B": tauc_AB, "C": tauc_C, "D": tauc_D}


def set_constants(config):
    # Sliding law set up so that tauc == beta**2:
    config.set_flag("basal_resistance.pseudo_plastic.enabled", True)
    config.set_number("basal_resistance.pseudo_plastic.q", 1.0)
    config.set_number("basal_resistance.pseudo_plastic.u_threshold",
                      convert(1.0, "m / s", "m / year"))

    # Isothermal Glen flow law (ISMIP-HOM uses a fixed ice softness).
    config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")
    config.set_number("stress_balance.blatter.Glen_exponent", 3.0)
    config.set_number("flow_law.isothermal_Glen.ice_softness",
                      convert(1e-16, "Pa-3 year-1", "Pa-3 s-1"))

    config.set_number("constants.ice.density", 910.0)
    config.set_number("constants.standard_gravity", 9.81)
    config.set_string("output.format", "netcdf4_serial")

    # Multigrid-compatible vertical grid: Mz = A*C^(N-1)+1. With coarsening
    # factor 2 and 3 mg levels (-bp_pc_mg_levels 3) this gives 9 -> 5 -> 3.
    config.set_number("stress_balance.blatter.Mz", 9)
    config.set_number("stress_balance.blatter.coarsening_factor", 2)


def init(testname, L):
    # the domain is a square [0, L] * [0, L]
    P = PISM.GridParameters(config)
    P.Lx = L / 2.0
    P.Ly = L / 2.0
    P.x0 = L / 2.0
    P.y0 = L / 2.0
    P.periodicity = PISM.XY_PERIODIC
    P.z = PISM.DoubleVector(np.linspace(0, 2000, 201))
    P.Mx = MX
    P.My = MY

    # use 3 grid points in the y direction for the x-z tests
    if testname in "BD":
        P.My = 3

    P.ownership_ranges_from_options(config, ctx.size)

    grid = PISM.Grid(ctx.ctx, P)

    geometry = PISM.Geometry(grid)
    grid.variables().add(geometry.ice_thickness)
    grid.variables().add(geometry.cell_type)
    grid.variables().add(geometry.bed_elevation)

    # enthalpy values are irrelevant: these tests use an isothermal flow law
    enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
    enthalpy.metadata(0).long_name("enthalpy of ice").units("J kg^-1")

    yield_stress = PISM.Scalar(grid, "tauc")
    yield_stress.metadata(0).long_name("basal yield stress").units("Pa")

    with PISM.vec.Access([yield_stress, geometry.ice_thickness, geometry.bed_elevation]):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            # Convert from "Pa year / m" to "Pa s / m"
            yield_stress[i, j] = tauc[testname](x, y, L) * seconds_per_year
            geometry.ice_thickness[i, j] = surface[testname](x, y, L) - bed[testname](x, y, L)
            geometry.bed_elevation[i, j] = bed[testname](x, y, L)

    geometry.sea_level_elevation.copy_from(geometry.bed_elevation)
    geometry.sea_level_elevation.shift(-100.0)   # ensure the whole domain is grounded
    geometry.ensure_consistency(0.0)

    return grid, geometry, enthalpy, yield_stress


def main():
    set_constants(config)
    grid, geometry, enthalpy, yield_stress = init(TEST, L)

    Mz = int(config.get_number("stress_balance.blatter.Mz"))
    cf = int(config.get_number("stress_balance.blatter.coarsening_factor"))

    model = PISM.BlatterISMIPHOM(grid, Mz, cf, tests[TEST])
    stress_balance = PISM.StressBalance(grid, model, PISM.BlatterMod(model))
    stress_balance.init()

    inputs = PISM.StressBalanceInputs()
    inputs.geometry = geometry
    inputs.basal_yield_stress = yield_stress
    inputs.enthalpy = enthalpy

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()

    # Time the solve itself (barriers so the timing is not skewed by load imbalance).
    comm.barrier()
    t0 = time.time()
    stress_balance.update(inputs, True)
    comm.barrier()
    solve_seconds = time.time() - t0

    if rank == 0:
        print("BENCH test=%s L_km=%.0f Mx=%d My=%d Mz=%d solve_seconds=%.2f" %
              (TEST, L / 1e3, grid.Mx(), grid.My(), Mz, solve_seconds))

    # Write the vertically-averaged velocity (ubar/vbar) -- a compact field that
    # is identical across solvers when they converge to the same answer -- as a
    # self-contained file for the correctness comparison.
    stress_balance.advective_velocity().dump(OUT)


if __name__ == "__main__":
    main()
