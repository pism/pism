#!/usr/bin/env python
"""This script runs the Blatter stress balance solver.
"""

import numpy as np
import PISM
import PISM.testing

ctx = PISM.Context()
config = ctx.config

def inputs_from_file(filename):
    grid = PISM.IceGrid.FromFile(ctx.ctx, filename, ["enthalpy", "thk"], PISM.CELL_CENTER)

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

    P = PISM.GridParameters(config)
    P.Lx = 10e3
    P.Ly = 10e3
    P.x0 = 0
    P.y0 = 0
    P.horizontal_size_from_options()
    P.vertical_grid_from_options(config)
    P.ownership_ranges_from_options(ctx.size)

    grid = PISM.IceGrid(ctx.ctx, P)

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

    yield_stress.set(10 * config.get_number("basal_yield_stress.constant.value"))

    return grid, geometry, enthalpy, yield_stress

if __name__ == "__main__":

    input_file = config.get_string("input.file")

    if input_file == "":
        grid, geometry, enthalpy, yield_stress = inputs_from_formulas()
    else:
        grid, geometry, enthalpy, yield_stress = inputs_from_file(input_file)

    Mz = int(config.get_number("stress_balance.blatter.Mz"))
    n_levels = 0

    model = PISM.Blatter(grid, Mz, n_levels)
    mod = PISM.BlatterMod(model)

    stress_balance = PISM.StressBalance(grid, model, mod)

    stress_balance.init()

    inputs = PISM.StressBalanceInputs()

    inputs.geometry = geometry
    inputs.basal_yield_stress = yield_stress
    inputs.enthalpy = enthalpy

    stress_balance.update(inputs, True)

    output_file = config.get_string("output.file_name")

    PISM.util.prepare_output(output_file)

    model.velocity_u_sigma().write(output_file)
    model.velocity_v_sigma().write(output_file)
    stress_balance.velocity_u().write(output_file)
    stress_balance.velocity_v().write(output_file)
    stress_balance.velocity_w().write(output_file)

    stress_balance.advective_velocity().write(output_file)

    stress_balance.basal_frictional_heating().write(output_file)
    stress_balance.volumetric_strain_heating().write(output_file)

    v_mag = PISM.IceModelVec2S(grid, "vel_mag", PISM.WITHOUT_GHOSTS)
    v_mag.set_attrs("diagnostic",
                    "magnitude of vertically-averaged horizontal velocity of ice",
                    "m second-1", "m year-1", "", 0)
    v_mag.set_to_magnitude(stress_balance.advective_velocity())

    v_mag.write(output_file)

    geometry.ice_thickness.write(output_file)
    geometry.bed_elevation.write(output_file)
    enthalpy.write(output_file)
    yield_stress.write(output_file)

    print(stress_balance.max_diffusivity())
