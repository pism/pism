#!/usr/bin/env python
"""This script reads input data from a file (set with "-i"), runs BlatterStressBalance,
and then saves results to a different file (set with "-o").

An input file has to contain

- ice thickness,
- bed elevation,
- basal yield stress,
- ice enthalpy

For now we assume that the sea level is zero.

"""

import PISM

ctx = PISM.Context()

config = ctx.config

input_file = config.get_string("input.file")

grid = PISM.IceGrid.FromFile(ctx.ctx, input_file, ["enthalpy"], PISM.CELL_CENTER)

geometry = PISM.Geometry(grid)

enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS)
enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

tauc     = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
tauc.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

geometry.ice_thickness.regrid(input_file, critical=True)
geometry.bed_elevation.regrid(input_file, critical=True)
geometry.sea_level_elevation.set(0.0)

geometry.ensure_consistency(0.0)

enthalpy.regrid(input_file, critical=True)
tauc.regrid(input_file, critical=True)

stress_balance = PISM.BlatterStressBalance(grid, ctx.enthalpy_converter)

stress_balance.init()

inputs = PISM.StressBalanceInputs()

inputs.geometry = geometry
inputs.basal_yield_stress = tauc
inputs.enthalpy = enthalpy

stress_balance.update(inputs, True)

output_file = config.get_string("output.file_name")

PISM.util.prepare_output(output_file)

stress_balance.velocity_u().write(output_file)
stress_balance.velocity_v().write(output_file)

stress_balance.velocity_u_sigma().write(output_file)
stress_balance.velocity_v_sigma().write(output_file)
