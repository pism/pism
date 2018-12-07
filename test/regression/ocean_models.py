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

config.set_string("ocean.delta_sl_2d.file", "delta_SL_input.nc")

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

# change the default melange back pressure fraction from 0 to 1. The default of zero makes
# it hard to test the modifier that scales this value.
config.set_double("ocean.constant.melange_back_pressure_fraction", 1.0)

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


def check_model(model, T, SMB, MBP):
    check(model.shelf_base_temperature(), T)
    check(model.shelf_base_mass_flux(), SMB)
    check(model.melange_back_pressure_fraction(), MBP)


def check_modifier(model, modifier, dT, dSMB, dMBP):
    check_difference(modifier.shelf_base_temperature(),
                     model.shelf_base_temperature(),
                     dT)

    check_difference(modifier.shelf_base_mass_flux(),
                     model.shelf_base_mass_flux(),
                     dSMB)

    check_difference(modifier.melange_back_pressure_fraction(),
                     model.melange_back_pressure_fraction(),
                     dMBP)


def constant_test():
    "Model Constant"

    depth = 1000.0                  # meters

    # compute mass flux
    melt_rate = config.get_double("ocean.constant.melt_rate", "m second-1")
    ice_density = config.get_double("constants.ice.density")
    mass_flux = melt_rate * ice_density

    # compute pressure melting temperature
    T0 = config.get_double("constants.fresh_water.melting_point_temperature")
    beta_CC = config.get_double("constants.ice.beta_Clausius_Clapeyron")
    g = config.get_double("constants.standard_gravity")

    pressure = ice_density * g * depth
    T_melting = T0 - beta_CC * pressure

    melange_back_pressure = 1.0

    grid = dummy_grid()
    geometry = create_geometry(grid)
    geometry.ice_thickness.set(depth)

    model = PISM.OceanConstant(grid)

    model.init(geometry)
    model.update(geometry, 0, 1)

    check_model(model, T_melting, mass_flux, melange_back_pressure)

    assert model.max_timestep(0).infinite() == True


def pik_test():
    "Model PIK"
    grid = dummy_grid()
    geometry = create_geometry(grid)

    depth = 1000.0                  # meters

    # compute pressure melting temperature
    ice_density = config.get_double("constants.ice.density")
    T0 = config.get_double("constants.fresh_water.melting_point_temperature")
    beta_CC = config.get_double("constants.ice.beta_Clausius_Clapeyron")
    g = config.get_double("constants.standard_gravity")

    pressure = ice_density * g * depth
    T_melting = T0 - beta_CC * pressure

    melange_back_pressure = 0.0

    mass_flux = 5.36591610659e-06  # stored mass flux value returned by the model

    # create the model
    geometry.ice_thickness.set(depth)

    model = PISM.OceanPIK(grid)

    model.init(geometry)
    model.update(geometry, 0, 1)

    check_model(model, T_melting, mass_flux, melange_back_pressure)

    assert model.max_timestep(0).infinite() == True


class GivenTest(TestCase):
    "Test the Given class"

    def setUp(self):
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.filename = "given_input.nc"

        self.temperature = 263.0
        self.mass_flux = 3e-3
        self.melange_back_pressure = 0.0

        create_given_input_file(self.filename, self.grid, self.temperature, self.mass_flux)

        config.set_string("ocean.given.file", self.filename)

    def runTest(self):
        "Model Given"

        model = PISM.OceanGiven(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        assert model.max_timestep(0).infinite() == True

        check_model(model, self.temperature, self.mass_flux, self.melange_back_pressure)

    def tearDown(self):
        os.remove(self.filename)


class GivenTHTest(TestCase):
    def setUp(self):

        depth = 1000.0
        salinity = 35.0
        potential_temperature = 270.0
        self.melange_back_pressure = 0.0
        self.temperature = 270.17909999999995
        self.mass_flux = -6.489250000000001e-05

        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.geometry.ice_thickness.set(depth)

        filename = "given_th_input.nc"
        self.filename = filename

        PISM.util.prepare_output(filename)

        Th = PISM.IceModelVec2S(self.grid, "theta_ocean", PISM.WITHOUT_GHOSTS)
        Th.set_attrs("climate", "potential temperature", "Kelvin", "")
        Th.set(potential_temperature)
        Th.write(filename)

        S = PISM.IceModelVec2S(self.grid, "salinity_ocean", PISM.WITHOUT_GHOSTS)
        S.set_attrs("climate", "ocean salinity", "g/kg", "")
        S.set(salinity)
        S.write(filename)

        config.set_string("ocean.th.file", self.filename)

    def runTest(self):
        "Model GivenTH"

        model = PISM.OceanGivenTH(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        assert model.max_timestep(0).infinite() == True

        check_model(model, self.temperature, self.mass_flux, self.melange_back_pressure)

    def tearDown(self):
        os.remove(self.filename)


class DeltaT(TestCase):
    def setUp(self):
        self.filename = "delta_T_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.OceanConstant(self.grid)
        self.dT = -5.0

        create_dummy_forcing_file(self.filename, "delta_T", "Kelvin", self.dT)

    def runTest(self):
        "Modifier Delta_T"

        modifier = PISM.OceanDeltaT(self.grid, self.model)

        options.setValue("-ocean_delta_T_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, self.dT, 0.0, 0.0)

    def tearDown(self):
        os.remove(self.filename)


class DeltaSMB(TestCase):
    def setUp(self):
        self.filename = "delta_SMB_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.OceanConstant(self.grid)
        self.dSMB = -5.0

        create_dummy_forcing_file(self.filename, "delta_mass_flux", "kg m-2 s-1", self.dSMB)

    def runTest(self):
        "Modifier Delta_SMB"

        modifier = PISM.OceanDeltaSMB(self.grid, self.model)

        options.setValue("-ocean_delta_mass_flux_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, 0.0, self.dSMB, 0.0)

    def tearDown(self):
        os.remove(self.filename)

class AnomalyBMB(TestCase):
    def setUp(self):
        self.filename = "delta_BMB_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.OceanConstant(self.grid)
        self.dBMB = -5.0

        delta_BMB = PISM.IceModelVec2S(self.grid, "shelf_base_mass_flux_anomaly",
                                       PISM.WITHOUT_GHOSTS)
        delta_BMB.set_attrs("climate_forcing",
                            "2D shelf base mass flux anomaly", "kg m-2 s-1", "")
        delta_BMB.set(self.dBMB)

        delta_BMB.dump(self.filename)

    def runTest(self):
        "Modifier Anomaly"

        config.set_string("ocean.anomaly.file", self.filename)

        modifier = PISM.OceanAnomaly(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, 0.0, self.dBMB, 0.0)

    def tearDown(self):
        os.remove(self.filename)

class FracMBP(TestCase):
    def setUp(self):
        self.filename = "frac_MBP_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.OceanConstant(self.grid)
        self.dMBP = 0.5

        create_dummy_forcing_file(self.filename, "frac_MBP", "1", self.dMBP)

    def runTest(self):
        "Modifier Frac_MBP"

        modifier = PISM.OceanFracMBP(self.grid, self.model)

        options.setValue("-ocean_frac_MBP_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        model = self.model

        check_difference(modifier.shelf_base_temperature(),
                         model.shelf_base_temperature(),
                         0.0)

        check_difference(modifier.shelf_base_mass_flux(),
                         model.shelf_base_mass_flux(),
                         0.0)

        check_ratio(modifier.melange_back_pressure_fraction(),
                    model.melange_back_pressure_fraction(),
                    self.dMBP)

    def tearDown(self):
        os.remove(self.filename)


class FracSMB(TestCase):
    def setUp(self):
        self.filename = "frac_SMB_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.OceanConstant(self.grid)
        self.dSMB = 0.5

        create_dummy_forcing_file(self.filename, "frac_mass_flux", "1", self.dSMB)

    def runTest(self):
        "Modifier Frac_SMB"

        modifier = PISM.OceanFracSMB(self.grid, self.model)

        options.setValue("-ocean_frac_mass_flux_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        model = self.model

        check_difference(modifier.shelf_base_temperature(),
                         model.shelf_base_temperature(),
                         0.0)

        check_ratio(modifier.shelf_base_mass_flux(),
                    model.shelf_base_mass_flux(),
                    self.dSMB)

        check_difference(modifier.melange_back_pressure_fraction(),
                         model.melange_back_pressure_fraction(),
                         0.0)

    def tearDown(self):
        os.remove(self.filename)


def create_delta_T_file(filename):
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


class Cache(TestCase):
    def setUp(self):
        self.filename = "dT.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.constant = PISM.OceanConstant(self.grid)
        self.delta_T = PISM.OceanDeltaT(self.grid, self.constant)

        create_delta_T_file(self.filename)
        options.setValue("-ocean_delta_T_file", self.filename)

        config.set_double("ocean.cache.update_interval", 2.0)

    def runTest(self):
        "Modifier Cache"

        modifier = PISM.OceanCache(self.grid, self.delta_T)

        modifier.init(self.geometry)

        t = 0
        dt = seconds_per_year

        diff = []
        while t < 4 * seconds_per_year:
            modifier.update(self.geometry, t, dt)

            original = sample(self.constant.shelf_base_temperature())
            cached = sample(modifier.shelf_base_temperature())

            diff.append(cached - original)

            t += dt

        np.testing.assert_almost_equal(diff, [1, 1, 3, 3])

    def tearDown(self):
        os.remove(self.filename)


class DeltaSL(TestCase):
    def setUp(self):
        self.filename = "delta_SL_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.SeaLevel(self.grid)
        self.dSL = -5.0

        create_dummy_forcing_file(self.filename, "delta_SL", "meters", self.dSL)

    def runTest(self):
        "Modifier Delta_SL"

        modifier = PISM.SeaLevelDelta(self.grid, self.model)

        options.setValue("-ocean_delta_sl_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_difference(modifier.elevation(),
                         self.model.elevation(),
                         self.dSL)

    def tearDown(self):
        os.remove(self.filename)


def create_delta_SL_file(filename, grid, sea_level_offset):
    PISM.util.prepare_output(filename)

    SL = PISM.IceModelVec2S(grid, "delta_SL", PISM.WITHOUT_GHOSTS)
    SL.set_attrs("forcing", "sea level forcing", "meters", "")
    SL.set(sea_level_offset)
    SL.write(filename)


class DeltaSL2D(TestCase):
    def setUp(self):
        self.filename = config.get_string("ocean.delta_sl_2d.file")
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.SeaLevel(self.grid)
        self.dSL = -5.0

        create_delta_SL_file(self.filename, self.grid, self.dSL)

    def runTest(self):
        "Modifier Delta_SL_2D"

        modifier = PISM.SeaLevelDelta2D(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_difference(modifier.elevation(),
                         self.model.elevation(),
                         self.dSL)

    def tearDown(self):
        os.remove(self.filename)
