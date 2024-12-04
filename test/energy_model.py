#!/usr/bin/env python3

import PISM
from PISM.util import convert

ctx = PISM.Context()

ctx.log.set_threshold(1)

def create_dummy_grid():
    "Create a dummy grid"
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.Grid(ctx.ctx, params)


def setup():
    global grid

    # "inputs" does not own its elements, so we need to make sure
    # these don't expire when setup() returns
    global basal_melt_rate, ice_thickness, inputs, surface_temp
    global climatic_mass_balance, basal_heat_flux, shelf_base_temp, cell_type
    global u, v, w, strain_heating3

    grid = create_dummy_grid()

    zero = PISM.Scalar(grid, "zero")
    zero.set(0.0)

    cell_type = PISM.CellType(grid, "mask")
    cell_type.set(PISM.MASK_GROUNDED)

    basal_heat_flux = PISM.Scalar(grid, "bheatflx")
    basal_heat_flux.set(convert(10, "mW m-2", "W m-2"))

    ice_thickness = PISM.model.createIceThicknessVec(grid)
    ice_thickness.set(4000.0)
    # TemperatureModel needs ice_thickness to set enthalpy in restart(...)
    grid.variables().add(ice_thickness)

    shelf_base_temp = PISM.Scalar(grid, "shelfbtemp")
    shelf_base_temp.set(260.0)

    surface_temp = PISM.Scalar(grid, "surface_temp")
    surface_temp.set(260.0)

    strain_heating3 = PISM.Array3D(grid, "sigma", PISM.WITHOUT_GHOSTS)

    u = PISM.Array3D(grid, "u", PISM.WITHOUT_GHOSTS)

    v = PISM.Array3D(grid, "v", PISM.WITHOUT_GHOSTS)

    w = PISM.Array3D(grid, "w", PISM.WITHOUT_GHOSTS)

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
    inputs.volumetric_heating_rate = strain_heating3
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
    except TypeError:
        pass

    print(model.stdout_flags())
    stats = model.stats()
    enthalpy = model.enthalpy()
    bmr = model.basal_melt_rate()


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

        print("Testing %s..." % M)

        print("* Bootstrapping using provided basal melt rate...")
        initialize(model)

        use_model(model)

        try:
            temperature = model.get_temperature()
        except:
            pass

        try:
            F = PISM.util.prepare_output("energy_model_state.nc")

            print("* Saving the model state...")
            model.write_model_state(F)

            print("* Restarting from a saved model state...")
            model.restart(F, 0)

            print("* Bootstrapping from a saved model state...")
            model.bootstrap(F,
                            ice_thickness,
                            surface_temp,
                            climatic_mass_balance,
                            basal_heat_flux)
        finally:
            F.close()


def test_temp_restart_from_enth():
    enth_model = PISM.EnthalpyModel(grid, None)
    temp_model = PISM.TemperatureModel(grid, None)

    initialize(enth_model)

    F = PISM.util.prepare_output("enth_model_state.nc")
    enth_model.write_model_state(F)

    temp_model.restart(F, 0)


def test_enth_restart_from_temp():
    enth_model = PISM.EnthalpyModel(grid, None)
    temp_model = PISM.TemperatureModel(grid, None)

    initialize(temp_model)

    F = PISM.util.prepare_output("temp_model_state.nc")
    temp_model.write_model_state(F)

    enth_model.restart(F, 0)


setup()

test_interface()
test_temp_restart_from_enth()
test_enth_restart_from_temp()
