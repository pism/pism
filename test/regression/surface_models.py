#!/usr/bin/env python
"""
Tests of PISM's surface models and modifiers.
"""

import PISM
from PISM.testing import *
import os
import numpy as np
from unittest import TestCase, SkipTest

config = PISM.Context().config

# reduce the grid size to speed this up
config.set_number("grid.Mx", 3)
config.set_number("grid.My", 5)
config.set_number("grid.Mz", 5)

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

# silence models' initialization messages
PISM.Context().log.set_threshold(1)

options = PISM.PETSc.Options()

def climatic_mass_balance(grid, value):
    SMB = PISM.IceModelVec2S(grid, "climatic_mass_balance", PISM.WITHOUT_GHOSTS)
    SMB.set_attrs("climate", "surface mass balance", "kg m-2 s-1",
                  "land_ice_surface_specific_mass_balance_flux")
    SMB.set(value)
    return SMB

def ice_surface_temp(grid, value):
    temperature = PISM.IceModelVec2S(grid, "ice_surface_temp", PISM.WITHOUT_GHOSTS)
    temperature.set_attrs("climate", "ice temperature at the top surface", "Kelvin", "")
    temperature.set(value)
    return temperature

def check_model(model, T, omega, SMB, mass=0.0, thickness=0.0):
    check(model.mass_flux(), SMB)
    check(model.temperature(), T)
    check(model.liquid_water_fraction(), omega)
    check(model.layer_mass(), mass)
    check(model.layer_thickness(), thickness)

def surface_simple(grid):
    return PISM.SurfaceSimple(grid, PISM.AtmosphereUniform(grid))

def check_modifier(model, modifier, T=0.0, omega=0.0, SMB=0.0, mass=0.0, thickness=0.0):
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

class DeltaT(TestCase):
    def setUp(self):
        self.filename = "delta_T_input.nc"
        self.grid = shallow_grid()
        self.model = surface_simple(self.grid)
        self.dT = -5.0
        self.geometry = PISM.Geometry(self.grid)

        create_scalar_forcing(self.filename, "delta_T", "Kelvin", [self.dT], [0])

    def test_surface_delta_t(self):
        "Modifier 'delta_T'"

        modifier = PISM.SurfaceDeltaT(self.grid, self.model)

        config.set_string("surface.delta_T.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT)

    def tearDown(self):
        os.remove(self.filename)

class LapseRates(TestCase):
    def setUp(self):
        self.filename = "reference_surface.nc"
        self.grid = shallow_grid()
        self.model = surface_simple(self.grid)
        self.dTdz = 1.0         # 1 Kelvin per km
        self.dT = -1.0
        self.dz = 1000.0

        self.geometry = PISM.Geometry(self.grid)

        # save current surface elevation to use it as a "reference" surface elevation
        self.geometry.ice_surface_elevation.dump(self.filename)

        config.set_string("surface.lapse_rate.file", self.filename)
        config.set_number("surface.lapse_rate.temperature_lapse_rate", self.dTdz)

    def test_surface_lapse_rate(self):
        "Modifier 'lapse_rate'"

        modifier = PISM.SurfaceLapseRates(self.grid, self.model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check that the temperature changed accordingly
        modifier.update(self.geometry, 0, 1)
        check_modifier(self.model, modifier, T=self.dT)

    def tearDown(self):
        os.remove(self.filename)

class Given(TestCase):
    def setUp(self):
        self.filename = "surface_given_input.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.SMB = 10.0
        self.T = 250.0

        output = PISM.util.prepare_output(self.filename)
        climatic_mass_balance(self.grid, self.SMB).write(output)
        ice_surface_temp(self.grid, self.T).write(output)

        config.set_string("surface.given.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_surface_given(self):
        "Model 'given'"

        atmosphere = PISM.AtmosphereUniform(self.grid)
        model = PISM.SurfaceFactory(self.grid, atmosphere).create("given")

        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        check_model(model, T=self.T, omega=0, SMB=self.SMB, mass=0.0, thickness=0.0)

def test_surface_elevation():
    "Model 'elevation'"
    T_min = -5.0
    T_max = 0.0
    z_min = 1000.0
    z_ela = 1100.0
    z_max = 1500.0
    M_min = -1.0
    M_max = 5.0

    grid = shallow_grid()
    geometry = PISM.Geometry(grid)
    geometry.ice_thickness.set(0.5 * (z_min + z_max))
    geometry.ensure_consistency(0.0)

    options.setValue("-ice_surface_temp", "{},{},{},{}".format(T_min, T_max, z_min, z_max))
    options.setValue("-climatic_mass_balance",
                     "{},{},{},{},{}".format(M_min, M_max, z_min, z_ela, z_max))

    T = PISM.util.convert(0.5 * (T_min + T_max), "Celsius", "Kelvin")
    SMB = PISM.util.convert(1.87504, "m/year", "m/s") * config.get_number("constants.ice.density")

    atmosphere = PISM.AtmosphereUniform(grid)
    model = PISM.SurfaceElevation(grid, atmosphere)

    model.init(geometry)

    model.update(geometry, 0, 1)

    check_model(model, T=T, SMB=SMB, omega=0, mass=0, thickness=0)

class TemperatureIndex(TestCase):
    def setUp(self):
        self.air_temp = config.get_number("atmosphere.uniform.temperature")
        self.precip = config.get_number("atmosphere.uniform.precipitation")

        self.grid = shallow_grid()

        self.geometry = PISM.Geometry(self.grid)
        # make sure that there's ice to melt
        self.geometry.ice_thickness.set(1000.0)

        T_above_zero = 1
        dt_days = 5
        self.T = 273.15 + T_above_zero
        self.dt = dt_days * 86400

        ice_density = config.get_number("constants.ice.density")
        beta_ice = config.get_number("surface.pdd.factor_ice")
        refreeze_fraction = config.get_number("surface.pdd.refreeze")
        PDD = dt_days * T_above_zero
        ice_melted = PDD * beta_ice
        refreeze = ice_melted * refreeze_fraction
        self.SMB = -(ice_melted - refreeze) * ice_density / self.dt

        config.set_number("atmosphere.uniform.temperature", self.T)
        # disable daily variability so that we can compute the number of PDDs exactly
        config.set_number("surface.pdd.std_dev", 0.0)
        # no precipitation
        config.set_number("atmosphere.uniform.precipitation", 0)

    def tearDown(self):
        config.set_number("atmosphere.uniform.temperature", self.air_temp)
        config.set_number("atmosphere.uniform.precipitation", self.precip)

    def test_surface_pdd(self):
        "Model 'pdd'"
        model = PISM.SurfaceTemperatureIndex(self.grid, PISM.AtmosphereUniform(self.grid))

        model.init(self.geometry)

        model.update(self.geometry, 0, self.dt)

        check_model(model, T=self.T, SMB=self.SMB, omega=0.0, mass=0.0, thickness=0.0)

class PIK(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.filename = "pik_input.nc"
        self.SMB = 10.0

        self.geometry.latitude.set(-80.0)
        self.geometry.ice_thickness.set(2000.0)
        self.geometry.ensure_consistency(0.0)

        output = PISM.util.prepare_output(self.filename)
        climatic_mass_balance(self.grid, self.SMB).write(output)

        config.set_string("input.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)
        config.set_string("input.file", "")

    def test_surface_pik(self):
        "Model 'pik'"

        model = PISM.SurfacePIK(self.grid, None)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, T=233.13, SMB=self.SMB, omega=0.0, mass=0.0, thickness=0.0)

def test_surface_simple():
    "Model 'simple'"
    grid = shallow_grid()
    geometry = PISM.Geometry(grid)

    atmosphere = PISM.AtmosphereUniform(grid)
    model = PISM.SurfaceSimple(grid, atmosphere)

    T = atmosphere.mean_annual_temp().numpy()[0,0]
    SMB = atmosphere.mean_precipitation().numpy()[0, 0]

    check_model(model, T=T, SMB=SMB, omega=0.0, mass=0.0, thickness=0.0)

class Anomaly(TestCase):
    def setUp(self):
        self.filename = "anomaly_input.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = surface_simple(self.grid)
        self.dT = -5.0
        self.dSMB = 20.0

        dT = PISM.IceModelVec2S(self.grid, "ice_surface_temp_anomaly", PISM.WITHOUT_GHOSTS)
        dT.set_attrs("climate", "temperature anomaly", "Kelvin", "")
        dT.set(self.dT)

        dSMB = PISM.IceModelVec2S(self.grid, "climatic_mass_balance_anomaly", PISM.WITHOUT_GHOSTS)
        dSMB.set_attrs("climate", "SMB anomaly", "kg m-2 s-1", "")
        dSMB.set(self.dSMB)

        output = PISM.util.prepare_output(self.filename)
        dT.write(output)
        dSMB.write(output)

        config.set_string("surface.anomaly.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_anomaly(self):
        "Modifier 'anomaly'"

        modifier = PISM.SurfaceAnomaly(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT, SMB=self.dSMB)

class Cache(TestCase):
    def setUp(self):
        self.filename = "dT.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.simple = surface_simple(self.grid)
        self.delta_T = PISM.SurfaceDeltaT(self.grid, self.simple)

        time_bounds = np.array([0, 1, 1, 2, 2, 3, 3, 4]) * seconds_per_year
        create_scalar_forcing(self.filename, "delta_T", "Kelvin", [1, 2, 3, 4],
                              times=None, time_bounds=time_bounds)

        config.set_string("surface.delta_T.file", self.filename)

        config.set_number("surface.cache.update_interval", 2.0)

    def test_surface_cache(self):
        "Modifier 'cache'"

        modifier = PISM.SurfaceCache(self.grid, self.delta_T)

        modifier.init(self.geometry)

        dt = seconds_per_year

        N = 4
        ts = np.arange(float(N)) * dt
        diff = []
        for t in ts:
            modifier.update(self.geometry, t, dt)

            original = sample(self.simple.temperature())
            cached = sample(modifier.temperature())

            diff.append(cached - original)

        np.testing.assert_almost_equal(diff, [1, 1, 3, 3])

    def tearDown(self):
        os.remove(self.filename)

class Forcing(TestCase):
    def setUp(self):
        self.filename = "forcing_input.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.simple = surface_simple(self.grid)

        config.set_string("surface.force_to_thickness_file", self.filename)

        ice_density = config.get_number("constants.ice.density")
        alpha = config.get_number("surface.force_to_thickness.alpha", "second-1")

        self.H = 500.0
        self.H_target = 1000.0
        self.dSMB = ice_density * alpha * (self.H_target - self.H)

        output = PISM.util.prepare_output(self.filename)
        # target thickness
        self.geometry.ice_thickness.set(self.H_target)
        self.geometry.ice_thickness.write(output)

        # ftt mask
        ftt_mask = PISM.IceModelVec2Int(self.grid, "ftt_mask", PISM.WITHOUT_GHOSTS)
        ftt_mask.set(1.0)
        ftt_mask.write(output)

        self.geometry.ice_thickness.set(self.H)
        self.geometry.ensure_consistency(0.0)

    def tearDown(self):
        os.remove(self.filename)

    def test_surface_forcing(self):
        "Modifier 'forcing'"

        modifier = PISM.SurfaceForceThickness(self.grid, self.simple)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_modifier(self.simple, modifier, SMB=self.dSMB)

class EISMINTII(TestCase):
    def test_eismint2(self):
        raise SkipTest("not implemented")

class Initialization(TestCase):
    def test_initialization(self):
        raise SkipTest("not implemented")
