#!/usr/bin/env python
"""This script runs the Blatter stress balance solver.
"""

import numpy as np
import PISM
import PISM.testing

ctx = PISM.Context()
config = ctx.config

def inputs_from_file(filename):
    grid = PISM.IceGrid.FromFile(ctx.ctx, filename, ["thk"], PISM.CELL_CENTER)

    geometry = PISM.Geometry(grid)

    geometry.ice_thickness.regrid(filename, critical=True)
    geometry.bed_elevation.regrid(filename, critical=True)
    geometry.sea_level_elevation.set(0.0)

    geometry.ensure_consistency(0.0)

    yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
    yield_stress.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

    yield_stress.regrid(filename, critical=False,
                        default_value=config.get_number("basal_yield_stress.constant.value"))

    enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS)
    enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

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

    return grid, geometry, enthalpy, yield_stress

def inputs_from_formulas():
    H_max = 1000.0

    grid = PISM.testing.shallow_grid(Mx=int(config.get_number("grid.Mx")),
                                     My=int(config.get_number("grid.My")),
                                     Lx=10e3,
                                     Ly=10e3)

    geometry = PISM.Geometry(grid)

    enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS)
    enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

    yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
    yield_stress.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

    with PISM.vec.Access(nocomm=[geometry.bed_elevation, geometry.ice_thickness]):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            geometry.bed_elevation[i, j] = 0.0
            geometry.ice_thickness[i, j] = np.sqrt(max(H_max**2 - (r / 10)**2, 0.0))

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

    return grid, geometry, enthalpy, yield_stress

if __name__ == "__main__":

    input_file = config.get_string("input.file")

    if input_file == "":
        grid, geometry, enthalpy, yield_stress = inputs_from_formulas()
    else:
        grid, geometry, enthalpy, yield_stress = inputs_from_file(input_file)

    Mz = int(config.get_number("stress_balance.blatter.Mz"))
    n_levels = 1

    stress_balance = PISM.Blatter(grid, Mz, n_levels)

    inputs = PISM.StressBalanceInputs()

    inputs.geometry = geometry
    inputs.basal_yield_stress = yield_stress
    inputs.enthalpy = enthalpy

    stress_balance.update(inputs, True)

    stress_balance.update(inputs, False)

    output_file = config.get_string("output.file_name")

    PISM.util.prepare_output(output_file)

    stress_balance.u_velocity().write(output_file)
    stress_balance.v_velocity().write(output_file)

    geometry.ice_thickness.write(output_file)
    geometry.bed_elevation.write(output_file)
    enthalpy.write(output_file)
    yield_stress.write(output_file)
