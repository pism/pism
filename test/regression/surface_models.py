#!/usr/bin/env python
"""
Tests of PISM's surface models and modifiers.
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


def create_geometry(grid):
    geometry = PISM.Geometry(grid)

    geometry.latitude.set(0.0)
    geometry.longitude.set(0.0)

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)

    geometry.ice_thickness.set(0.0)
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


def dummy_grid():
    "Create a dummy grid"
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)


def check(vec, value):
    "Check if values of vec are almost equal to value."
    np.testing.assert_almost_equal(sample(vec), value)


def check_difference(A, B, value):
    "Check if the difference between A and B is almost equal to value."
    np.testing.assert_almost_equal(sample(A) - sample(B), value)


def check_model(model, T, omega, SMB,
                mass=0.0, thickness=0.0, accumulation=0.0, melt=0.0, runoff=0.0):
    check(model.mass_flux(), SMB)
    check(model.temperature(), T)
    check(model.liquid_water_fraction(), omega)
    check(model.layer_mass(), mass)
    check(model.layer_thickness(), thickness)
    check(model.accumulation(), accumulation)
    check(model.melt(), melt)
    check(model.runoff(), runoff)


def check_modifier(model, modifier,
                   T=0.0, omega=0.0, SMB=0.0, mass=0.0, thickness=0.0,
                   accumulation=0.0, melt=0.0, runoff=0.0):
    check_difference(modifier.mass_flux(),
                     model.mass_flux(),
                     SMB)

    check_difference(modifier.temperature(),
                     model.temperature(),
                     T)

    check_difference(modifier.liquid_water_fraction(),
                     model.liquid_water_fraction(),
                     omega)

    check_difference(modifier.layer_mass(),
                     model.layer_mass(),
                     mass)

    check_difference(modifier.layer_thickness(),
                     model.layer_thickness(),
                     thickness)

    check_difference(modifier.accumulation(),
                     model.accumulation(),
                     accumulation)

    check_difference(modifier.melt(),
                     model.melt(),
                     melt)

    check_difference(modifier.runoff(),
                     model.runoff(),
                     runoff)

class Given(TestCase):
    "Model given"

    def create_given_input_file(self, filename, grid, temperature, mass_flux):
        PISM.util.prepare_output(filename)

        T = PISM.IceModelVec2S(grid, "ice_surface_temp", PISM.WITHOUT_GHOSTS)
        T.set_attrs("climate", "ice surface temperature", "Kelvin", "")
        T.set(temperature)
        T.write(filename)

        M = PISM.IceModelVec2S(grid, "climatic_mass_balance", PISM.WITHOUT_GHOSTS)
        M.set_attrs("climate", "top surface mass balance", "kg m-2 s-1", "")
        M.set(mass_flux)
        M.write(filename)

    def setUp(self):
        self.filename = "given_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.T = 272.15
        self.M = 1001.0

        self.create_given_input_file(self.filename, self.grid, self.T, self.M)

    def runTest(self):
        atmosphere = PISM.AtmosphereUniform(self.grid)

        config.set_string("surface.given.file", self.filename)

        model = PISM.SurfaceGiven(self.grid, atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, self.T, 0.0, self.M, accumulation=self.M)

    def tearDown(self):
        os.remove(self.filename)

class DeltaT(TestCase):
    def setUp(self):
        self.filename = "delta_T_input.nc"
        self.grid = dummy_grid()
        self.model = PISM.SurfaceEISMINTII(self.grid, ord('A'))
        self.dT = -5.0
        self.geometry = create_geometry(self.grid)

        create_dummy_forcing_file(self.filename, "delta_T", "Kelvin", self.dT)

    def runTest(self):
        "Modifier Delta_T"

        modifier = PISM.SurfaceDeltaT(self.grid, self.model)

        options.setValue("-surface_delta_T_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT)

    def tearDown(self):
        os.remove(self.filename)


class LapseRates(TestCase):
    def setUp(self):
        self.filename = "reference_surface.nc"
        self.grid = dummy_grid()
        self.model = PISM.SurfaceEISMINTII(self.grid, ord('A'))
        self.dTdz = 1.0         # 1 Kelvin per km
        self.dT = -1.0
        self.dz = 1000.0

        self.geometry = create_geometry(self.grid)

        # save current surface elevation to use it as a "reference" surface elevation
        self.geometry.ice_surface_elevation.dump(self.filename)

        config.set_string("surface.lapse_rate.file", self.filename)

        options.setValue("-temp_lapse_rate", self.dTdz)

    def runTest(self):
        "Modifier lapse_rate"

        modifier = PISM.SurfaceLapseRates(self.grid, self.model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check that the temperature changed accordingly
        modifier.update(self.geometry, 0, 1)
        check_modifier(self.model, modifier, T=self.dT)

    def tearDown(self):
        os.remove(self.filename)

class Elevation(TestCase):
    def setUp(self):
        pass
    def runTest(self):
        raise NotImplementedError
    def tearDown(self):
        pass

class TemperatureIndex(TestCase):
    def setUp(self):
        pass
    def runTest(self):
        raise NotImplementedError
    def tearDown(self):
        pass

class PIK(TestCase):
    def create_input(self, filename, grid, mass_flux):
        PISM.util.prepare_output(filename)

        M = PISM.IceModelVec2S(grid, "climatic_mass_balance", PISM.WITHOUT_GHOSTS)
        M.set_attrs("climate", "top surface mass balance", "kg m-2 s-1", "")
        M.set(mass_flux)
        M.write(filename)

    def setUp(self):
        self.filename = "pik_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.M = 1001.0
        self.T = 273.15 + 30.0

        self.create_input(self.filename, self.grid, self.M)

        config.set_string("input.file", self.filename)

    def runTest(self):
        model = PISM.SurfacePIK(self.grid, PISM.AtmosphereUniform(self.grid))

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, self.T, 0.0, self.M, accumulation=self.M)

    def tearDown(self):
        os.remove(self.filename)

class Simple(TestCase):
    def setUp(self):
        self.grid = dummy_grid()
        self.atmosphere = PISM.AtmosphereUniform(self.grid)
        self.geometry = create_geometry(self.grid)

    def runTest(self):
        atmosphere = self.atmosphere

        model = PISM.SurfaceSimple(self.grid, atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        T = sample(atmosphere.mean_annual_temp())
        M = sample(atmosphere.mean_precipitation())

        check_model(model, T, 0.0, M, accumulation=M)

    def tearDown(self):
        pass

class Anomaly(TestCase):
    def setUp(self):
        pass
    def runTest(self):
        raise NotImplementedError
    def tearDown(self):
        pass

class Cache(TestCase):
    def setUp(self):
        pass
    def runTest(self):
        raise NotImplementedError
    def tearDown(self):
        pass

class ForceThickness(TestCase):
    def setUp(self):
        pass
    def runTest(self):
        raise NotImplementedError
    def tearDown(self):
        pass

if __name__ == "__main__":

    t = Given()

    t.setUp()
    t.runTest()
    t.tearDown()
