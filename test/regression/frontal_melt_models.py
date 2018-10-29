#!/usr/bin/env python
"""
Tests of PISM's ocean models and modifiers.
"""

import PISM
import sys
import os
import numpy as np
from unittest import TestCase
import netCDF4

config = PISM.Context().config

# reduce the grid size to speed this up
config.set_double("grid.Mx", 3)
config.set_double("grid.My", 3)
config.set_double("grid.Mz", 5)

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

log = PISM.Context().log
# silence models' initialization messages
log.set_threshold(1)

options = PISM.PETSc.Options()


def frontal_melt_from_discharge_and_thermal_forcing(h, Qsg, TF):

    A = config.get_double("frontal_melt.parameter_a")
    B = config.get_double("frontal_melt.parameter_b")
    alpha = config.get_double("frontal_melt.power_alpha")
    beta = config.get_double("frontal_melt.power_beta")

    return (A * h * Qsg ** alpha + B) * TF ** beta


def create_geometry(grid):
    geometry = PISM.Geometry(grid)

    geometry.cell_area.set(grid.dx() * grid.dy())
    geometry.latitude.set(0.0)
    geometry.longitude.set(0.0)

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)

    geometry.ice_thickness.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)

    geometry.ensure_consistency(0.0)

    return geometry


def sample(vec):
    with PISM.vec.Access(nocomm=[vec]):
        return vec[0, 0]


def create_dummy_forcing_file(filename, variable_name, units, value):
    f = netCDF4.Dataset(filename, "w")
    f.createDimension("time", 1)
    t = f.createVariable("time", "d", ("time",))
    t.units = "seconds"
    delta_T = f.createVariable(variable_name, "d", ("time",))
    delta_T.units = units
    t[0] = 0.0
    delta_T[0] = value
    f.close()


def dummy_grid():
    "Create a dummy grid"
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)


def create_given_input_file(filename, grid, temperature, mass_flux):
    PISM.util.prepare_output(filename)

    T = PISM.IceModelVec2S(grid, "shelfbtemp", PISM.WITHOUT_GHOSTS)
    T.set_attrs("climate", "shelf base temperature", "Kelvin", "")
    T.set(temperature)
    T.write(filename)

    M = PISM.IceModelVec2S(grid, "shelfbmassflux", PISM.WITHOUT_GHOSTS)
    M.set_attrs("climate", "shelf base mass flux", "kg m-2 s-1", "")
    M.set(mass_flux)
    M.write(filename)


def check(vec, value):
    "Check if values of vec are almost equal to value."
    np.testing.assert_almost_equal(sample(vec), value)


def check_difference(A, B, value):
    "Check if the difference between A and B is almost equal to value."
    np.testing.assert_almost_equal(sample(A) - sample(B), value)


def check_ratio(A, B, value):
    "Check if the ratio of A and B is almost equal to value."
    b = sample(B)
    if b != 0:
        np.testing.assert_almost_equal(sample(A) / b, value)
    else:
        np.testing.assert_almost_equal(sample(A), 0.0)


def check_model(model, melt_rate):
    check(model.frontal_melt_rate(), melt_rate)


def check_modifier(model, modifier, dT, dSMB, dMBP):
    check_difference(modifier.shelf_base_temperature(), model.shelf_base_temperature(), dT)

    check_difference(modifier.shelf_base_mass_flux(), model.shelf_base_mass_flux(), dSMB)

    check_difference(modifier.melange_back_pressure_fraction(), model.melange_back_pressure_fraction(), dMBP)


def constant_test():
    "Model Constant"

    # compute mass flux
    melt_rate = config.get_double("frontal_melt.constant.melt_rate", "m second-1")

    grid = dummy_grid()
    geometry = create_geometry(grid)

    model = PISM.FrontalMeltConstant(grid)

    model.init(geometry)
    model.update(geometry, 0, 1)

    check_model(model, melt_rate)

    assert model.max_timestep(0).infinite() == True


class DischargeRoutingTest(TestCase):
    def setUp(self):

        depth = 1000.0
        self.potential_temperature = 4.0
        self.subglacial_discharge = 1.0  # check units

        self.grid = dummy_grid()

        Qsg = PISM.IceModelVec2S(self.grid, "tendency_of_subglacial_water_mass_at_grounding_line", PISM.WITHOUT_GHOSTS)
        Qsg.set_attrs("climate", "subglacial discharge at grounding line", "m day-1", "m s-1")
        Qsg.set(self.subglacial_discharge)

        # This seg faults but why?
        self.grid.variables().add(Qsg)

        self.geometry = create_geometry(self.grid)
        self.geometry.ice_thickness.set(depth)

        filename = "discharge_routing_input.nc"
        self.filename = filename

        PISM.util.prepare_output(filename)

        Th = PISM.IceModelVec2S(self.grid, "theta_ocean", PISM.WITHOUT_GHOSTS)
        Th.set_attrs("climate", "potential temperature", "Kelvin", "")
        Th.set(self.potential_temperature)
        Th.write(filename)

        config.set_string("frontal_melt.routing.file", self.filename)
        config.set_string("hydrology.model", "routing")
        melt_rate(self.depth, self.subglacial_discharge, self.potential_temperature)
        print(melt_rate)

    def runTest(self):
        "Model DischargeRouting"

        model = PISM.FrontalMeltDischargeRouting(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        assert model.max_timestep(0).infinite() == True

        check_model(model, melt_rate)

    def tearDown(self):
        os.remove(self.filename)


class GivenTest(TestCase):
    def setUp(self):

        self.frontal_melt_rate = 100.0

        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        filename = "given_input.nc"
        self.filename = filename

        PISM.util.prepare_output(filename)

        Fmr = PISM.IceModelVec2S(self.grid, "frontalmeltrate", PISM.WITHOUT_GHOSTS)
        Fmr.set_attrs("climate", "frontal melt rate", "m/year", "")
        Fmr.set(self.frontal_melt_rate)

        config.set_string("frontal_melt.given.file", self.filename)

    def runTest(self):
        "Model Given"

        model = PISM.FrontalMeltGiven(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        assert model.max_timestep(0).infinite() == True

        check_model(model, self.temperature, self.mass_flux, self.melange_back_pressure)

    def tearDown(self):
        os.remove(self.filename)
