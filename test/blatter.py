#!/usr/bin/env python
"""This script runs the Blatter stress balance solver.
"""

import PISM
import PISM.testing

ctx = PISM.Context()

config = ctx.config

grid = PISM.testing.shallow_grid(Mx=int(config.get_number("grid.Mx")),
                                 My=int(config.get_number("grid.My")),
                                 Lx=10e3,
                                 Ly=10e3)

geometry = PISM.Geometry(grid)

enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS)
enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

yield_stress     = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
yield_stress.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

with PISM.vec.Access(nocomm=geometry.bed_elevation):
    for (i, j) in grid.points():
        geometry.bed_elevation[i, j] = (grid.x(i) + grid.Lx()) * 1e-4

geometry.ice_thickness.set(1000.0)
geometry.sea_level_elevation.set(0.0)

geometry.ensure_consistency(0.0)

ice_surface_temp = PISM.IceModelVec2S(grid, "ice_surface_temp", PISM.WITHOUT_GHOSTS)
ice_surface_temp.set(250.0)

smb = PISM.IceModelVec2S(grid, "smb", PISM.WITHOUT_GHOSTS)
smb.set(0.0)

basal_heat_flux = PISM.IceModelVec2S(grid, "basal_heat_flux", PISM.WITHOUT_GHOSTS)
basal_heat_flux.set(0.0)

PISM.bootstrap_ice_enthalpy(geometry.ice_thickness,
                            ice_surface_temp,
                            smb,
                            basal_heat_flux,
                            enthalpy)

yield_stress.set(config.get_number("basal_yield_stress.constant.value"))

stress_balance = PISM.BlatterStressBalance(grid, ctx.enthalpy_converter)

stress_balance.init()

inputs = PISM.StressBalanceInputs()

inputs.geometry = geometry
inputs.basal_yield_stress = yield_stress
inputs.enthalpy = enthalpy

stress_balance.update(inputs, True)

output_file = config.get_string("output.file_name")

PISM.util.prepare_output(output_file)

stress_balance.velocity_u().write(output_file)
stress_balance.velocity_v().write(output_file)

stress_balance.velocity_u_sigma().write(output_file)
stress_balance.velocity_v_sigma().write(output_file)

geometry.ice_thickness.write(output_file)
geometry.bed_elevation.write(output_file)
enthalpy.write(output_file)
yield_stress.write(output_file)
