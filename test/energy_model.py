#!/usr/bin/env python

import PISM
from PISM.util import convert

ctx = PISM.Context()


def create_dummy_grid():
    "Create a dummy grid"
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)


def setup():
    global grid

    # "inputs" does not own its elements, so we need to make sure
    # these don't expire when setup() returns
    global basal_melt_rate, ice_thickness, inputs, surface_temp
    global climatic_mass_balance, basal_heat_flux, shelf_base_temp, cell_type
    global u, v, w, strain_heating3

    grid = create_dummy_grid()

    zero = PISM.IceModelVec2S()
    zero.create(grid, "zero", PISM.WITHOUT_GHOSTS)
    zero.set(0.0)

    cell_type = PISM.IceModelVec2CellType()
    cell_type.create(grid, "mask", PISM.WITHOUT_GHOSTS)
    cell_type.set(PISM.MASK_GROUNDED)

    basal_heat_flux = PISM.IceModelVec2S()
    basal_heat_flux.create(grid, "bheatflx", PISM.WITHOUT_GHOSTS)
    basal_heat_flux.set(convert(10, "mW m-2", "W m-2"))

    ice_thickness = PISM.model.createIceThicknessVec(grid)
    ice_thickness.set(4000.0)
    # TemperatureModel needs ice_thickness to set enthalpy in restart(...)
    grid.variables().add(ice_thickness)

    shelf_base_temp = PISM.IceModelVec2S()
    shelf_base_temp.create(grid, "shelfbtemp", PISM.WITHOUT_GHOSTS)
    shelf_base_temp.set(260.0)

    surface_temp = PISM.IceModelVec2S()
    surface_temp.create(grid, "surface_temp", PISM.WITHOUT_GHOSTS)
    surface_temp.set(260.0)

    strain_heating3 = PISM.IceModelVec3()
    strain_heating3.create(grid, "sigma", PISM.WITHOUT_GHOSTS)

    u = PISM.IceModelVec3()
    u.create(grid, "u", PISM.WITHOUT_GHOSTS)

    v = PISM.IceModelVec3()
    v.create(grid, "v", PISM.WITHOUT_GHOSTS)

    w = PISM.IceModelVec3()
    w.create(grid, "w", PISM.WITHOUT_GHOSTS)

    ice_thickness.set(4000.0)
    u.set(0.0)
    v.set(0.0)
    w.set(0.0)

    basal_melt_rate = zero
    climatic_mass_balance = zero

    inputs = PISM.EnergyModelInputs()

    inputs.cell_type = cell_type
    inputs.basal_frictional_heating = zero
    inputs.basal_heat_flux = basal_heat_flux
    inputs.ice_thickness = ice_thickness
    inputs.surface_liquid_fraction = zero
    inputs.shelf_base_temp = shelf_base_temp
    inputs.surface_temp = surface_temp
    inputs.till_water_thickness = zero
    inputs.strain_heating3 = strain_heating3
    inputs.u3 = u
    inputs.v3 = v
    inputs.w3 = w

dt = convert(1, "years", "seconds")

def use_model(model):
    print("* Performing a time step...")
    model.update(0, dt, inputs)

    try:
        model.update(0, dt)
        raise Exception("this should fail")
    except RuntimeError:
        pass

    try:
        model.max_timestep(0)
        raise Exception("this should fail")
    except RuntimeError:
        pass

    print(model.stdout_flags())
    stats = model.stats()
    enthalpy = model.get_enthalpy()
    bmr = model.get_basal_melt_rate()


def initialize(model):
    model.initialize(basal_melt_rate,
                     ice_thickness,
                     surface_temp,
                     climatic_mass_balance,
                     basal_heat_flux)


def test_interface():
    "Use the EnergyModel interface to make sure the code runs."
    for M in [PISM.EnthalpyModel, PISM.DummyEnergyModel, PISM.TemperatureModel]:
        model = M(grid, None)

        print("")
        print("Testing %s..." % M)

        print("* Bootstrapping using provided basal melt rate...")
        initialize(model)

        use_model(model)

        try:
            temperature = model.get_temperature()
        except:
            pass

        pio = PISM.util.prepare_output("energy_model_state.nc")

        print("* Saving the model state...")
        model.write_model_state(pio)

        print("* Restarting from a saved model state...")
        model.restart(pio, 0)

        print("* Bootstrapping from a saved model state...")
        model.bootstrap(pio,
                        ice_thickness,
                        surface_temp,
                        climatic_mass_balance,
                        basal_heat_flux)

        pio.close()


def test_temp_restart_from_enth():
    enth_model = PISM.EnthalpyModel(grid, None)
    temp_model = PISM.TemperatureModel(grid, None)

    initialize(enth_model)

    pio = PISM.util.prepare_output("enth_model_state.nc")
    enth_model.write_model_state(pio)

    temp_model.restart(pio, 0)


def test_enth_restart_from_temp():
    enth_model = PISM.EnthalpyModel(grid, None)
    temp_model = PISM.TemperatureModel(grid, None)

    initialize(temp_model)

    pio = PISM.util.prepare_output("temp_model_state.nc")
    temp_model.write_model_state(pio)

    enth_model.restart(pio, 0)


setup()

test_interface()
test_temp_restart_from_enth()
test_enth_restart_from_temp()
