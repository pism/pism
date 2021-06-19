#!/usr/bin/env python3
"""
Tests of PISM's ocean models and modifiers.
"""

import PISM
from PISM.testing import *
from PISM.testing import filename as tmp_name
import os
import numpy as np
from unittest import TestCase

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

def water_column_pressure(ice_thickness):
    "Average column pressure"
    ice_density   = config.get_number("constants.ice.density")
    water_density = config.get_number("constants.sea_water.density")
    water_depth   = ice_thickness * ice_density / water_density
    g = config.get_number("constants.standard_gravity")
    return 0.5 * water_density * g * water_depth**2 / ice_thickness

def create_given_input_file(filename, grid, temperature, mass_flux):
    PISM.util.prepare_output(filename)

    T = PISM.IceModelVec2S(grid, "shelfbtemp", PISM.WITHOUT_GHOSTS)
    T.set_attrs("climate", "shelf base temperature", "Kelvin", "Kelvin", "", 0)
    T.set(temperature)
    T.write(filename)

    M = PISM.IceModelVec2S(grid, "shelfbmassflux", PISM.WITHOUT_GHOSTS)
    M.set_attrs("climate", "shelf base mass flux", "kg m-2 s-1", "kg m-2 s-1", "", 0)
    M.set(mass_flux)
    M.write(filename)

def check_model(model, T, SMB, water_column_pressure):
    "Check values returned by an ocean model"
    check(model.shelf_base_temperature(), T)
    check(model.shelf_base_mass_flux(), SMB)
    check(model.average_water_column_pressure(), water_column_pressure)

def check_modifier(model, modifier, dT, dSMB, dwater_column_pressure):
    check_difference(modifier.shelf_base_temperature(),
                     model.shelf_base_temperature(),
                     dT)

    check_difference(modifier.shelf_base_mass_flux(),
                     model.shelf_base_mass_flux(),
                     dSMB)

    check_difference(modifier.average_water_column_pressure(),
                     model.average_water_column_pressure(),
                     dwater_column_pressure)

def constant_test():
    "Model Constant"

    ice_thickness = 1000.0                  # meters

    # compute mass flux
    melt_rate = config.get_number("ocean.constant.melt_rate", "m second-1")
    ice_density = config.get_number("constants.ice.density")
    mass_flux = melt_rate * ice_density

    # compute pressure melting temperature
    T0 = config.get_number("constants.fresh_water.melting_point_temperature")
    beta_CC = config.get_number("constants.ice.beta_Clausius_Clapeyron")
    g = config.get_number("constants.standard_gravity")

    pressure = ice_density * g * ice_thickness
    T_melting = T0 - beta_CC * pressure

    average_water_column_pressure = water_column_pressure(ice_thickness)

    grid = shallow_grid()
    geometry = PISM.Geometry(grid)
    geometry.ice_thickness.set(ice_thickness)
    geometry.bed_elevation.set(-2 * ice_thickness)
    geometry.ensure_consistency(0.0)

    model = PISM.OceanConstant(grid)

    model.init(geometry)
    model.update(geometry, 0, 1)

    check_model(model, T_melting, mass_flux, average_water_column_pressure)

    assert model.max_timestep(0).infinite() == True

def pik_test():
    "Model PIK"
    grid = shallow_grid()
    geometry = PISM.Geometry(grid)

    ice_thickness = 1000.0                  # meters

    # compute pressure melting temperature
    ice_density = config.get_number("constants.ice.density")
    T0 = config.get_number("constants.fresh_water.melting_point_temperature")
    beta_CC = config.get_number("constants.ice.beta_Clausius_Clapeyron")
    g = config.get_number("constants.standard_gravity")

    pressure = ice_density * g * ice_thickness
    T_melting = T0 - beta_CC * pressure

    average_water_column_pressure = water_column_pressure(ice_thickness)

    mass_flux = 5.36591610659e-06  # stored mass flux value returned by the model

    # create the model
    geometry.ice_thickness.set(ice_thickness)
    geometry.bed_elevation.set(-2 * ice_thickness)
    geometry.ensure_consistency(0.0)

    model = PISM.OceanPIK(grid)

    model.init(geometry)
    model.update(geometry, 0, 1)

    check_model(model, T_melting, mass_flux, average_water_column_pressure)

    assert model.max_timestep(0).infinite() == True

class GivenTest(TestCase):
    "Test the Given class"

    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.filename = tmp_name("ocean_given_input")

        self.temperature = 263.0
        self.mass_flux = 3e-3

        ice_thickness = 1000.0
        self.geometry.ice_thickness.set(ice_thickness)
        self.geometry.bed_elevation.set(-2 * ice_thickness)
        self.geometry.ensure_consistency(0.0)

        self.average_water_column_pressure = water_column_pressure(ice_thickness)

        create_given_input_file(self.filename, self.grid, self.temperature, self.mass_flux)

        config.set_string("ocean.given.file", self.filename)

    def test_ocean_given(self):
        "Model Given"

        model = PISM.OceanGiven(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        assert model.max_timestep(0).infinite() == True

        check_model(model, self.temperature, self.mass_flux,
                    self.average_water_column_pressure)

    def tearDown(self):
        os.remove(self.filename)

class GivenTHTest(TestCase):
    def setUp(self):

        ice_thickness = 1000.0
        salinity = 35.0
        potential_temperature = 270.0
        self.average_water_column_pressure = water_column_pressure(ice_thickness)
        self.temperature = 270.17909999999995
        self.mass_flux = -6.489250000000001e-05

        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.geometry.ice_thickness.set(ice_thickness)
        self.geometry.bed_elevation.set(-2 * ice_thickness)
        self.geometry.ensure_consistency(0.0)

        self.filename = tmp_name("ocean_given_th_input")
        filename = self.filename

        PISM.util.prepare_output(filename)

        Th = PISM.IceModelVec2S(self.grid, "theta_ocean", PISM.WITHOUT_GHOSTS)
        Th.set_attrs("climate", "potential temperature", "Kelvin", "Kelvin", "", 0)
        Th.set(potential_temperature)
        Th.write(filename)

        S = PISM.IceModelVec2S(self.grid, "salinity_ocean", PISM.WITHOUT_GHOSTS)
        S.set_attrs("climate", "ocean salinity", "g/kg", "g/kg", "", 0)
        S.set(salinity)
        S.write(filename)

        config.set_string("ocean.th.file", self.filename)

    def test_ocean_th(self):
        "Model GivenTH"

        model = PISM.OceanGivenTH(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        assert model.max_timestep(0).infinite() == True

        check_model(model, self.temperature, self.mass_flux, self.average_water_column_pressure)

    def tearDown(self):
        os.remove(self.filename)

class DeltaT(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_delta_T_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)

        self.model = PISM.OceanConstant(self.grid)
        self.dT = -5.0

        PISM.testing.create_scalar_forcing(self.filename, "delta_T", "Kelvin", [self.dT], [0])

    def test_ocean_delta_t(self):
        "Modifier Delta_T"

        modifier = PISM.OceanDeltaT(self.grid, self.model)

        config.set_string("ocean.delta_T.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, self.dT, 0.0, 0.0)

    def tearDown(self):
        os.remove(self.filename)

class DeltaSMB(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_delta_SMB_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)
        self.model = PISM.OceanConstant(self.grid)
        self.dSMB = -5.0

        create_scalar_forcing(self.filename, "delta_mass_flux", "kg m-2 s-1",
                              [self.dSMB], [0])

    def test_ocean_delta_smb(self):
        "Modifier Delta_SMB"

        modifier = PISM.OceanDeltaSMB(self.grid, self.model)

        config.set_string("ocean.delta_mass_flux.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, 0.0, self.dSMB, 0.0)

    def tearDown(self):
        os.remove(self.filename)

class AnomalyBMB(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_delta_BMB_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)

        self.model = PISM.OceanConstant(self.grid)
        self.dBMB = -5.0

        delta_BMB = PISM.IceModelVec2S(self.grid, "shelf_base_mass_flux_anomaly",
                                       PISM.WITHOUT_GHOSTS)
        delta_BMB.set_attrs("climate_forcing",
                            "2D shelf base mass flux anomaly", "kg m-2 s-1", "kg m-2 s-1",
                            "", 0)
        delta_BMB.set(self.dBMB)

        delta_BMB.dump(self.filename)

    def test_ocean_anomaly(self):
        "Modifier Anomaly"

        config.set_string("ocean.anomaly.file", self.filename)

        modifier = PISM.OceanAnomaly(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, 0.0, self.dBMB, 0.0)

    def tearDown(self):
        os.remove(self.filename)

class DeltaMBP(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_delta_MBP_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.model = PISM.OceanConstant(self.grid)
        self.dP = 100.0         # Pascal (over the thickness of melange)
        self.H = 1000.0         # meters
        self.dP_ice = config.get_number("ocean.delta_MBP.melange_thickness") * self.dP / self.H

        self.geometry.ice_thickness.set(self.H)
        self.geometry.bed_elevation.set(-2 * self.H)
        self.geometry.ensure_consistency(0.0)

        create_scalar_forcing(self.filename, "delta_MBP", "Pa", [self.dP], [0])

    def test_ocean_delta_mpb(self):
        "Modifier Delta_MBP"

        modifier = PISM.OceanDeltaMBP(self.grid, self.model)

        config.set_string("ocean.delta_MBP.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        model = self.model

        check_modifier(model, modifier, 0.0, 0.0, self.dP_ice)

    def tearDown(self):
        os.remove(self.filename)

class FracMBP(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_frac_MBP_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.model = PISM.OceanConstant(self.grid)
        self.Lambda = 0.5
        self.H = 1000.0

        g = config.get_number("constants.standard_gravity")
        rho_w = config.get_number("constants.sea_water.density")
        rho_i = config.get_number("constants.ice.density")

        def P(depth, density, g):
            return 0.5 * density * g * depth

        P_water = water_column_pressure(self.H)
        P_ice = 0.5 * rho_i * g * self.H
        self.dP = self.Lambda * (P_ice - P_water)

        self.geometry.ice_thickness.set(self.H)
        self.geometry.bed_elevation.set(-2 * self.H)
        self.geometry.ensure_consistency(0.0)

        create_scalar_forcing(self.filename, "frac_MBP", "1", [self.Lambda], [0])

    def test_ocean_frac_mpb(self):
        "Modifier Frac_MBP"

        modifier = PISM.OceanFracMBP(self.grid, self.model)

        config.set_string("ocean.frac_MBP.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        model = self.model

        check_difference(modifier.shelf_base_temperature(),
                         model.shelf_base_temperature(),
                         0.0)

        check_difference(modifier.shelf_base_mass_flux(),
                         model.shelf_base_mass_flux(),
                         0.0)

        check_difference(modifier.average_water_column_pressure(),
                         model.average_water_column_pressure(),
                         self.dP)

    def tearDown(self):
        os.remove(self.filename)


class FracSMB(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_frac_SMB_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)
        self.model = PISM.OceanConstant(self.grid)
        self.dSMB = 0.5

        create_scalar_forcing(self.filename, "frac_mass_flux", "1", [self.dSMB], [0])

    def test_ocean_frac_smb(self):
        "Modifier Frac_SMB"

        modifier = PISM.OceanFracSMB(self.grid, self.model)

        config.set_string("ocean.frac_mass_flux.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        model = self.model

        check_difference(modifier.shelf_base_temperature(),
                         model.shelf_base_temperature(),
                         0.0)

        check_ratio(modifier.shelf_base_mass_flux(),
                    model.shelf_base_mass_flux(),
                    self.dSMB)

        check_difference(modifier.average_water_column_pressure(),
                         model.average_water_column_pressure(),
                         0.0)

    def tearDown(self):
        os.remove(self.filename)

class Cache(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_dT")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)

        self.constant = PISM.OceanConstant(self.grid)
        self.delta_T = PISM.OceanDeltaT(self.grid, self.constant)

        time_bounds = np.array([0, 1, 1, 2, 2, 3, 3, 4]) * seconds_per_year
        create_scalar_forcing(self.filename, "delta_T", "Kelvin", [1, 2, 3, 4],
                              times=None, time_bounds=time_bounds)

        config.set_string("ocean.delta_T.file", self.filename)
        config.set_number("ocean.cache.update_interval", 2.0)

    def test_ocean_cache(self):
        "Modifier Cache"

        model = self.constant
        modifier = PISM.OceanCache(self.grid, self.delta_T)

        modifier.init(self.geometry)

        dt = seconds_per_year

        N = 4
        ts = np.arange(float(N)) * dt
        diff = []
        for t in ts:
            modifier.update(self.geometry, t, dt)

            check_difference(model.shelf_base_mass_flux(),
                             modifier.shelf_base_mass_flux(), 0.0)

            check_difference(model.average_water_column_pressure(),
                             modifier.average_water_column_pressure(), 0.0)

            original = sample(model.shelf_base_temperature())
            cached = sample(modifier.shelf_base_temperature())

            diff.append(cached - original)

        np.testing.assert_almost_equal(diff, [1, 1, 3, 3])

    def tearDown(self):
        os.remove(self.filename)

class DeltaSL(TestCase):
    def setUp(self):
        self.filename = tmp_name("ocean_delta_SL_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)
        self.model = PISM.SeaLevel(self.grid)
        self.dSL = -5.0

        create_scalar_forcing(self.filename, "delta_SL", "meters", [self.dSL], [0])

    def test_ocean_delta_sl(self):
        "Modifier Delta_SL"

        modifier = PISM.SeaLevelDelta(self.grid, self.model)

        config.set_string("ocean.delta_sl.file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_difference(modifier.elevation(),
                         self.model.elevation(),
                         self.dSL)

    def tearDown(self):
        os.remove(self.filename)

class DeltaSL2D(TestCase):
    def create_delta_SL_file(self, filename, times, sea_level_offsets):
        output = PISM.util.prepare_output(filename, append_time=False)

        SL = PISM.IceModelVec2S(self.grid, "delta_SL", PISM.WITHOUT_GHOSTS)
        SL.set_attrs("forcing", "sea level forcing", "meters", "meters", "", 0)

        for k in range(len(times)):
            PISM.append_time(output, config, times[k])
            SL.set(sea_level_offsets[k])
            SL.write(output)

    def setUp(self):
        self.filename = tmp_name("ocean_delta_SL_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)
        self.model = PISM.SeaLevel(self.grid)
        self.dSL = -5.0

        self.create_delta_SL_file(self.filename, [0, seconds_per_year], [0, 2 * self.dSL])

        config.set_string("ocean.delta_sl_2d.file", self.filename)

    def test_ocean_delta_sl_2d(self):
        "Modifier Delta_SL_2D"

        modifier = PISM.SeaLevelDelta2D(self.grid, self.model)

        modifier.init(self.geometry)
        # Use a one second time step to try to sample sea level forcing midway through the
        # interval from 0 to 1 year.
        modifier.update(self.geometry, 0.5 * seconds_per_year, 1)
        check_difference(modifier.elevation(),
                         self.model.elevation(),
                         self.dSL)

    def tearDown(self):
        os.remove(self.filename)


if __name__ == "__main__":
    PISM.Context().log.set_threshold(3)

    t = DeltaMBP()
    t.setUp()
    t.test_ocean_delta_mpb()
    t.tearDown()
