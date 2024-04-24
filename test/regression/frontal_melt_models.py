#!/usr/bin/env python3
"""
Tests of PISM's frontal melt models.
"""

import PISM
from PISM.util import convert
import sys, os, numpy
from unittest import TestCase
import netCDF4

config = PISM.Context().config

# reduce the grid size to speed this up
config.set_number("grid.Mx", 3)
config.set_number("grid.My", 3)
config.set_number("grid.Mz", 5)

log = PISM.Context().log
# silence models' initialization messages
log.set_threshold(1)

options = PISM.PETSc.Options()

seconds_per_day = 86400

def create_geometry(grid):
    geometry = PISM.Geometry(grid)

    geometry.latitude.set(0.0)
    geometry.longitude.set(0.0)

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)

    geometry.ice_thickness.set(1000.0)
    geometry.ice_area_specific_volume.set(0.0)

    geometry.ensure_consistency(0.0)

    return geometry

def sample(vec):
    return vec.numpy()[0,0]

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

def create_grid():
    "Create a dummy grid"
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.Grid(ctx.ctx, params)

def create_given_input_file(filename, grid, temperature, mass_flux):
    PISM.util.prepare_output(filename)

    T = PISM.Scalar(grid, "shelfbtemp")
    T.set_attrs("climate", "shelf base temperature", "kelvin", "kelvin", "", 0)
    T.set(temperature)
    T.write(filename)

    M = PISM.Scalar(grid, "shelfbmassflux")
    M.set_attrs("climate", "shelf base mass flux", "kg m-2 s-1", "kg m-2 s-1", "", 0)
    M.set(mass_flux)
    M.write(filename)

def check(vec, value):
    "Check if values of vec are almost equal to value."
    numpy.testing.assert_almost_equal(sample(vec), value)

def check_model(model, melt_rate):
    check(model.frontal_melt_rate(), melt_rate)

def constant_test():
    "Model Constant"

    # compute mass flux
    melt_rate = config.get_number("frontal_melt.constant.melt_rate", "m second-1")

    grid = create_grid()
    geometry = create_geometry(grid)

    inputs = PISM.FrontalMeltInputs()
    water_flux = PISM.Scalar(grid, "water_flux")
    water_flux.set(0.0)
    inputs.geometry = geometry
    inputs.subglacial_water_flux = water_flux

    model = PISM.FrontalMeltConstant(grid)

    model.init(geometry)
    model.update(inputs, 0, 1)

    check_model(model, melt_rate)

    assert model.max_timestep(0).infinite()

class DischargeRoutingTest(TestCase):

    def frontal_melt(self, h, q_sg, TF):
        """
        h:    water depth, meters
        q_sg: subglacial water flux, m / day
        TF:   thermal forcing, Celsius

        Returns the melt rate in m / day
        """
        alpha = config.get_number("frontal_melt.routing.power_alpha")
        beta  = config.get_number("frontal_melt.routing.power_beta")
        A     = config.get_number("frontal_melt.routing.parameter_a")
        B     = config.get_number("frontal_melt.routing.parameter_b")

        return (A * h * q_sg ** alpha + B) * TF ** beta

    def setUp(self):
        self.depth = 1000.0              # meters
        self.potential_temperature = 4.0 # Celsius
        self.water_flux = 10.0           # m / day

        self.grid = create_grid()

        self.theta = PISM.Scalar(self.grid, "theta_ocean")
        self.theta.set(self.potential_temperature)

        self.Qsg = PISM.Scalar(self.grid, "subglacial_water_flux")
        self.Qsg.metadata(0).long_name("subglacial water flux").units("m2 / s")

        grid_spacing = 0.5 * (self.grid.dx() + self.grid.dy())
        cross_section_area = self.depth * grid_spacing

        self.Qsg.set(self.water_flux * cross_section_area / (grid_spacing * seconds_per_day))

        self.geometry = create_geometry(self.grid)
        self.geometry.ice_thickness.set(self.depth)
        self.geometry.sea_level_elevation.set(self.depth)

        self.geometry.ensure_consistency(config.get_number("geometry.ice_free_thickness_standard"))

        self.inputs = PISM.FrontalMeltInputs()
        self.inputs.geometry = self.geometry
        self.inputs.subglacial_water_flux = self.Qsg

    def runTest(self):
        "Model DischargeRouting"

        model = PISM.FrontalMeltDischargeRouting(self.grid)

        model.initialize(self.theta)

        model.update(self.inputs, 0, 1)

        melt_rate = self.frontal_melt(self.depth, self.water_flux, self.potential_temperature)
        # convert from m / day to m / s
        melt_rate /= seconds_per_day

        check_model(model, melt_rate)

        assert model.max_timestep(0).infinite()

    def tearDown(self):
        pass

class GivenTest(TestCase):
    def create_input(self, filename, melt_rate):
        PISM.util.prepare_output(filename)

        Fmr = PISM.Scalar(self.grid, "frontal_melt_rate")
        Fmr.metadata(0).long_name("frontal melt rate").units("m / s")

        Fmr.set(melt_rate)

        Fmr.write(filename)

    def setUp(self):

        self.frontal_melt_rate = 100.0

        self.grid = create_grid()
        self.geometry = create_geometry(self.grid)

        self.filename = "given_input.nc"

        self.create_input(self.filename, self.frontal_melt_rate)

        config.set_string("frontal_melt.given.file", self.filename)

        self.water_flux = PISM.Scalar(self.grid, "water_flux")
        self.water_flux.set(0.0)

        self.inputs = PISM.FrontalMeltInputs()
        self.inputs.geometry = self.geometry
        self.inputs.subglacial_water_flux = self.water_flux

    def runTest(self):
        "Model Given"

        model = PISM.FrontalMeltGiven(self.grid)
        model.init(self.geometry)

        model.update(self.inputs, 0, 1)

        assert model.max_timestep(0).infinite()

        check_model(model, self.frontal_melt_rate)

    def tearDown(self):
        os.remove(self.filename)

if __name__ == "__main__":

    t = DischargeRoutingTest()

    t.setUp()
    t.runTest()
    t.tearDown()
