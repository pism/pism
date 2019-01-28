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
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.atmosphere = PISM.AtmosphereUniform(self.grid)

    def runTest(self):

        config.set_string("surface.pdd.method", "expectation_integral")

        model = PISM.SurfaceTemperatureIndex(self.grid, self.atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        T = config.get_double("atmosphere.uniform.temperature", "Kelvin")
        omega = 0.0

        accumulation = config.get_double("atmosphere.uniform.precipitation", "kg m-2 second-1")
        melt         = accumulation
        runoff       = melt * (1.0 - config.get_double("surface.pdd.refreeze"))
        SMB          = accumulation - runoff

        check_model(model, T, omega, SMB, accumulation=accumulation, melt=melt,
                    runoff=runoff)

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
        config.set_string("input.file", "")

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
        self.filename = "surface_anomaly_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))
        self.dSMB = -5.0
        self.dT = 2.0

        PISM.util.prepare_output(self.filename)

        delta_SMB = PISM.IceModelVec2S(self.grid, "climatic_mass_balance_anomaly",
                                       PISM.WITHOUT_GHOSTS)
        delta_SMB.set_attrs("climate_forcing",
                            "2D surface mass flux anomaly", "kg m-2 s-1", "")
        delta_SMB.set(self.dSMB)

        delta_SMB.write(self.filename)

        delta_T = PISM.IceModelVec2S(self.grid, "ice_surface_temp_anomaly",
                                     PISM.WITHOUT_GHOSTS)
        delta_T.set_attrs("climate_forcing",
                          "2D surface temperature anomaly", "Kelvin", "")
        delta_T.set(self.dT)

        delta_T.write(self.filename)

    def runTest(self):
        "Modifier Anomaly"

        config.set_string("surface.anomaly.file", self.filename)

        modifier = PISM.SurfaceAnomaly(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT, SMB=self.dSMB)

    def tearDown(self):
        os.remove(self.filename)

class Cache(TestCase):
    def create_delta_T_file(self, filename):
        """Create a delta_T input file covering the interval [0, 4] years. When used with the
        'Cache' modifier with the update interval of 2, this data will be sampled every 2
        years, producing [1, 1, 3, 3].
        """
        f = netCDF4.Dataset(filename, "w")
        f.createDimension("time", 4)
        f.createDimension("bnds", 2)
        t = f.createVariable("time", "d", ("time",))
        t.units = "seconds"
        t.bounds = "time_bounds"
        t_bnds = f.createVariable("time_bounds", "d", ("time", "bnds"))
        delta_T = f.createVariable("delta_T", "d", ("time",))
        delta_T.units = "Kelvin"
        t[:] = [0.5, 1.5, 2.5, 3.5]
        t[:] *= seconds_per_year
        t_bnds[:, 0] = [0, 1, 2, 3]
        t_bnds[:, 1] = [1, 2, 3, 4]
        t_bnds[:, :] *= seconds_per_year
        delta_T[:] = [1, 2, 3, 4]
        f.close()

    def setUp(self):
        self.filename = "dT.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.constant = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))
        self.delta_T = PISM.SurfaceDeltaT(self.grid, self.constant)

        self.create_delta_T_file(self.filename)
        options.setValue("-surface_delta_T_file", self.filename)

        config.set_double("surface.cache.update_interval", 2.0)

    def runTest(self):
        "Modifier Cache"

        modifier = PISM.SurfaceCache(self.grid, self.delta_T)

        modifier.init(self.geometry)

        t = 0
        dt = seconds_per_year

        diff = []
        while t < 4 * seconds_per_year:
            modifier.update(self.geometry, t, dt)

            original = sample(self.constant.temperature())
            cached = sample(modifier.temperature())

            diff.append(cached - original)

            t += dt

        np.testing.assert_almost_equal(diff, [1, 1, 3, 3])

    def tearDown(self):
        os.remove(self.filename)

class ForceThickness(TestCase):
    def setUp(self):
        pass
    def runTest(self):
        raise NotImplementedError
    def tearDown(self):
        pass

if __name__ == "__main__":

    t = TemperatureIndex()

    t.setUp()
    t.runTest()
    t.tearDown()
