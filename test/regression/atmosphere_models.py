#!/usr/bin/env python3
"""
Tests of PISM's atmosphere models and modifiers.
"""

import PISM
from PISM.testing import *
from PISM.testing import filename as tmp_name
import os
import numpy as np
from unittest import TestCase

# set run duration to 1 second so that all forcing used here spans the duration of the run
time = PISM.Context().time
time.set_start(0)
time.set_end(1)

config = PISM.Context().config

# reduce the grid size to speed this up
config.set_number("grid.Mx", 3)
config.set_number("grid.My", 5) # non-square grid
config.set_number("grid.Mz", 2)

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

# silence models' initialization messages
PISM.Context().log.set_threshold(1)

def write_state(model):
    "Test writing of the model state"

    o_filename = tmp_name("atmosphere_model_state")
    o_diagnostics = tmp_name("atmosphere_diagnostics")

    try:
        output = PISM.util.prepare_output(o_filename)
        model.define_model_state(output)
        model.write_model_state(output)
        output.close()

        ds = model.diagnostics()
        output = PISM.util.prepare_output(o_diagnostics)

        for d in ds:
            ds[d].define(output, PISM.PISM_DOUBLE)

        for d in ds:
            ds[d].compute().write(output)

        output.close()

    finally:
        os.remove(o_filename)
        os.remove(o_diagnostics)

def check_model(model, T, P, ts=None, Ts=None, Ps=None):
    check(model.air_temperature(), T)
    check(model.precipitation(), P)

    model.init_timeseries(ts)

    try:
        model.begin_pointwise_access()
        np.testing.assert_almost_equal(model.temp_time_series(0, 0), Ts)
        np.testing.assert_almost_equal(model.precip_time_series(0, 0), Ps)
    finally:
        model.end_pointwise_access()

    write_state(model)

    model.max_timestep(ts[0])

def check_modifier(model, modifier, T=0.0, P=0.0, ts=None, Ts=None, Ps=None):
    check_difference(modifier.air_temperature(),
                     model.air_temperature(),
                     T)
    check_difference(modifier.precipitation(),
                     model.precipitation(),
                     P)

    model.init_timeseries(ts)
    modifier.init_timeseries(ts)

    try:
        model.begin_pointwise_access()
        modifier.begin_pointwise_access()

        Ts_model = np.array(model.temp_time_series(0, 0))
        Ts_modifier = np.array(modifier.temp_time_series(0, 0))

        Ps_model = np.array(model.precip_time_series(0, 0))
        Ps_modifier = np.array(modifier.precip_time_series(0, 0))

        np.testing.assert_almost_equal(Ts_modifier - Ts_model, Ts)
        np.testing.assert_almost_equal(Ps_modifier - Ps_model, Ps)
    finally:
        modifier.end_pointwise_access()
        model.end_pointwise_access()

    write_state(modifier)

    modifier.max_timestep(ts[0])

def precipitation(grid, value):
    precip = PISM.IceModelVec2S(grid, "precipitation")
    precip.set_attrs("climate", "precipitation", "kg m-2 s-1", "kg m-2 s-1",
                     "precipitation_flux", 0)
    precip.set(value)
    return precip

def air_temperature(grid, value):
    temperature = PISM.IceModelVec2S(grid, "air_temp")
    temperature.set_attrs("climate", "near-surface air temperature", "Kelvin", "Kelvin", "", 0)
    temperature.set(value)
    return temperature

class PIK(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_pik_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.geometry.latitude.set(-80.0)

        self.P = 10.0           # this is very high, but that's fine

        precipitation(self.grid, self.P).dump(self.filename)

        config.set_string("atmosphere.pik.file", self.filename)

    def test_atmosphere_pik(self):
        "Model 'pik'"

        parameterizations = {"martin" : (248.13, [248.13], [self.P]),
                             "huybrechts_dewolde" : (252.59, [237.973373], [self.P]),
                             "martin_huybrechts_dewolde" : (248.13, [233.51337298], [self.P]),
                             "era_interim" : (256.27, [243.0939774], [self.P]),
                             "era_interim_sin" : (255.31577, [241.7975841], [self.P]),
                             "era_interim_lon" : (248.886139, [233.3678998], [self.P])}

        for p, (T, Ts, Ps) in parameterizations.items():
            config.set_string("atmosphere.pik.parameterization", p)

            model = PISM.AtmospherePIK(self.grid)
            model.init(self.geometry)

            # t and dt are irrelevant here
            model.update(self.geometry, 0, 1)

            check_model(model, T=T, P=self.P, ts=[0.5], Ts=Ts, Ps=Ps)

        assert model.max_timestep(0).infinite()

        try:
            config.set_string("atmosphere.pik.parameterization", "invalid")
            model = PISM.AtmospherePIK(self.grid)
            assert False, "failed to catch an invalid parameterization"
        except RuntimeError:
            pass

    def tearDown(self):
        os.remove(self.filename)

class DeltaT1D(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_delta_T_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dT = -5.0

        create_scalar_forcing(self.filename, "delta_T", "Kelvin",
                              [self.dT], [0], time_bounds=[0, 1])

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_delta_t(self):
        "Modifier Delta_T"

        config.set_string("atmosphere.delta_T.file", self.filename)

        modifier = PISM.AtmosphereDeltaT(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT, ts=[0.5], Ts=[self.dT], Ps=[0])

class DeltaT2D(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_delta_T_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.delta_T = 5.0

        delta_T = PISM.IceModelVec2S(self.grid, "delta_T")
        delta_T.set_attrs("climate", "temperature offset", "K", "K", "", 0)
        delta_T.set(self.delta_T)

        try:
            output = PISM.util.prepare_output(self.filename)
            delta_T.write(output)
        finally:
            output.close()

        config.set_string("atmosphere.delta_T.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_frac_p(self):
        "Modifier 'delta_T': 2D offsets"

        modifier = PISM.AtmosphereDeltaT(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.delta_T, P=0.0,
                       ts=[0.5], Ts=[self.delta_T], Ps=[0])

class DeltaP1D(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_delta_P_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dP = 5.0

        create_scalar_forcing(self.filename, "delta_P", "kg m-2 s-1",
                              [self.dP], [0], time_bounds=[0, 1])

        config.set_string("atmosphere.delta_P.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_delta_p(self):
        "Modifier 'delta_P'"

        modifier = PISM.AtmosphereDeltaP(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, P=self.dP, ts=[0.5], Ts=[0], Ps=[self.dP])

class DeltaP2D(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_delta_P_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.delta_P = 5.0

        delta_P = PISM.IceModelVec2S(self.grid, "delta_P")
        delta_P.set_attrs("climate", "precipitation offset",
                          "kg m-2 s-1", "kg m-2 s-1", "", 0)
        delta_P.set(self.delta_P)

        try:
            output = PISM.util.prepare_output(self.filename)
            delta_P.write(output)
        finally:
            output.close()

        config.set_string("atmosphere.delta_P.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_frac_p(self):
        "Modifier 'delta_P': 2D offsets"

        modifier = PISM.AtmosphereDeltaP(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, P=self.delta_P, T=0.0,
                       ts=[0.5], Ps=[self.delta_P], Ts=[0])

class Given(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_given_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.P = 10.0
        self.T = 250.0

        output = PISM.util.prepare_output(self.filename)
        precipitation(self.grid, self.P).write(output)
        air_temperature(self.grid, self.T).write(output)

        config.set_string("atmosphere.given.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_given(self):
        "Model 'given'"

        model = PISM.AtmosphereFactory(self.grid).create("given")

        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        check_model(model, T=self.T, P=self.P, ts=[0.5], Ts=[self.T], Ps=[self.P])

class SeaRISE(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_searise_input")
        self.grid = shallow_grid()

        self.geometry = PISM.Geometry(self.grid)
        self.geometry.latitude.set(70.0)
        self.geometry.longitude.set(-45.0)
        self.geometry.ice_thickness.set(2500.0)
        self.geometry.ensure_consistency(0.0)

        self.P = 10.0

        output = PISM.util.prepare_output(self.filename)
        precipitation(self.grid, self.P).write(output)

        config.set_string("atmosphere.searise_greenland.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_searise_greenland(self):
        "Model 'searise_greenland'"

        model = PISM.AtmosphereSeaRISEGreenland(self.grid)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, P=self.P, T=251.9085, ts=[0.5], Ts=[238.66192632], Ps=[self.P])

class YearlyCycle(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_yearly_cycle")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.T_mean = 250.0
        self.T_summer = 270.0
        self.P = 15.0

        output = PISM.util.prepare_output(self.filename)
        precipitation(self.grid, self.P).write(output)

        T_mean = PISM.IceModelVec2S(self.grid, "air_temp_mean_annual")
        T_mean.set_attrs("climate", "mean annual near-surface air temperature", "K", "K", "", 0)
        T_mean.set(self.T_mean)
        T_mean.write(output)

        T_summer = PISM.IceModelVec2S(self.grid, "air_temp_mean_summer")
        T_summer.set_attrs("climate", "mean summer near-surface air temperature", "K", "K", "", 0)
        T_summer.set(self.T_summer)
        T_summer.write(output)

        config.set_string("atmosphere.yearly_cycle.file", self.filename)

        # FIXME: test "-atmosphere_yearly_cycle_scaling_file", too

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_yearly_cycle(self):
        "Model 'yearly_cycle'"

        model = PISM.AtmosphereCosineYearlyCycle(self.grid)

        model.init(self.geometry)

        one_year = 365 * 86400.0
        model.update(self.geometry, 0, one_year)

        summer_peak = config.get_number("atmosphere.fausto_air_temp.summer_peak_day") / 365.0

        ts = np.linspace(0, one_year, 13)
        cycle = np.cos(2.0 * np.pi * (ts / one_year - summer_peak))
        T = (self.T_summer - self.T_mean) * cycle + self.T_mean
        P = np.zeros_like(T) + self.P

        check_model(model, T=self.T_mean, P=self.P, ts=ts, Ts=T, Ps=P)

class OneStation(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_one_station")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.T = 263.15
        self.P = 10.0

        time_name = config.get_string("time.dimension_name")

        output = PISM.util.prepare_output(self.filename, append_time=True)

        output.redef()

        output.define_dimension("nv", 2)
        output.define_variable("time_bounds", PISM.PISM_DOUBLE, [time_name, "nv"])
        output.write_attribute(time_name, "bounds", "time_bounds")

        output.define_variable("precipitation", PISM.PISM_DOUBLE, [time_name])
        output.write_attribute("precipitation", "units", "kg m-2 s-1")

        output.define_variable("air_temp", PISM.PISM_DOUBLE, [time_name])
        output.write_attribute("air_temp", "units", "Kelvin")

        output.write_variable("precipitation", [0], [1], [self.P])
        output.write_variable("air_temp", [0], [1], [self.T])
        output.write_variable("time_bounds", [0, 0], [1, 2], [0, 1])

        output.close()

        config.set_string("atmosphere.one_station.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_one_station(self):
        "Model 'weather_station'"

        model = PISM.AtmosphereWeatherStation(self.grid)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, P=self.P, T=self.T, ts=[0.5], Ts=[self.T], Ps=[self.P])

class Uniform(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.P = 5.0
        self.T = 250.0

        config.set_number("atmosphere.uniform.temperature", self.T)
        config.set_number("atmosphere.uniform.precipitation", self.P)

    def test_atmosphere_uniform(self):
        "Model 'uniform'"
        model = PISM.AtmosphereUniform(self.grid)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        P = PISM.util.convert(self.P, "kg m-2 year-1", "kg m-2 s-1")
        check_model(model, T=self.T, P=P, ts=[0.5], Ts=[self.T], Ps=[P])

class Anomaly(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_anomaly_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dT = -5.0
        self.dP = 20.0

        dT = PISM.IceModelVec2S(self.grid, "air_temp_anomaly")
        dT.set_attrs("climate", "air temperature anomaly", "Kelvin", "Kelvin", "", 0)
        dT.set(self.dT)

        dP = PISM.IceModelVec2S(self.grid, "precipitation_anomaly")
        dP.set_attrs("climate", "precipitation anomaly", "kg m-2 s-1", "kg m-2 s-1", "", 0)
        dP.set(self.dP)

        output = PISM.util.prepare_output(self.filename)
        dT.write(output)
        dP.write(output)

        config.set_string("atmosphere.anomaly.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_anomaly(self):
        "Modifier 'anomaly'"

        modifier = PISM.AtmosphereAnomaly(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT, P=self.dP,
                       ts=[0.5], Ts=[self.dT], Ps=[self.dP])

class PrecipScaling(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_precip_scaling_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dT = 5.0

        create_scalar_forcing(self.filename, "delta_T", "Kelvin",
                              [self.dT], [0], time_bounds=[0, 1])

        config.set_string("atmosphere.precip_scaling.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_precip_scaling(self):
        "Modifier 'precip_scaling'"

        modifier = PISM.AtmospherePrecipScaling(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, P=1.3373514942327523e-05,
                       ts=[0.5], Ts=[0], Ps=[1.33735149e-05])

class FracP1D(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_frac_P_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.P_ratio = 5.0

        create_scalar_forcing(self.filename, "frac_P", "1",
                              [self.P_ratio], [0], time_bounds=[0, 1])

        config.set_string("atmosphere.frac_P.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_frac_p(self):
        "Modifier 'frac_P': 1D scaling"

        modifier = PISM.AtmosphereFracP(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_ratio(modifier.precipitation(), self.model.precipitation(),
                    self.P_ratio)

        check_modifier(self.model, modifier, T=0, P=0.00012675505856327396,
                       ts=[0.5], Ts=[0], Ps=[0.00012676])

class FracP2D(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_frac_P_input")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.P_ratio = 5.0

        frac_P = PISM.IceModelVec2S(self.grid, "frac_P")
        frac_P.set_attrs("climate", "precipitation scaling", "1", "1", "", 0)
        frac_P.set(self.P_ratio)

        try:
            output = PISM.util.prepare_output(self.filename)
            frac_P.write(output)
        finally:
            output.close()

        config.set_string("atmosphere.frac_P.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_frac_p(self):
        "Modifier 'frac_P': 2D scaling"

        modifier = PISM.AtmosphereFracP(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        check_ratio(modifier.precipitation(), self.model.precipitation(),
                    self.P_ratio)

        check_modifier(self.model, modifier, T=0, P=0.00012675505856327396,
                       ts=[0.5], Ts=[0], Ps=[0.00012676])


class ElevationChange(TestCase):
    def setUp(self):
        self.filename = tmp_name("atmosphere_reference_surface")
        self.grid = shallow_grid()
        self.dTdz = 1.0         # Kelvin per km
        self.dPdz = 1000.0      # (kg/m^2)/year per km
        self.dz = 1000.0        # m
        self.dT = -self.dTdz * self.dz / 1000.0
        self.dP = -PISM.util.convert(self.dPdz * self.dz / 1000.0, "kg m-2 year-1", "kg m-2 s-1")
        self.precip_dTdz = 2.0  # Kelvin per km

        self.geometry = PISM.Geometry(self.grid)

        # save current surface elevation to use it as a "reference" surface elevation
        self.geometry.ice_surface_elevation.dump(self.filename)

        config.set_string("atmosphere.elevation_change.file", self.filename)

        config.set_number("atmosphere.elevation_change.precipitation.lapse_rate", self.dPdz)

        config.set_number("atmosphere.elevation_change.precipitation.temp_lapse_rate", self.precip_dTdz)

        config.set_number("atmosphere.elevation_change.temperature_lapse_rate", self.dTdz)

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_elevation_change_shift(self):
        "Modifier 'elevation_change': lapse rate"

        config.set_string("atmosphere.elevation_change.precipitation.method", "shift")

        model = PISM.AtmosphereUniform(self.grid)
        modifier = PISM.AtmosphereElevationChange(self.grid, model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check that the temperature changed accordingly
        modifier.update(self.geometry, 0, 1)
        check_modifier(model, modifier, T=self.dT, P=self.dP,
                       ts=[0.5], Ts=[self.dT], Ps=[self.dP])

    def test_atmosphere_elevation_change_scale(self):
        "Modifier 'elevation_change': scaling"

        config.set_string("atmosphere.elevation_change.precipitation.method", "scale")

        model = PISM.AtmosphereUniform(self.grid)
        modifier = PISM.AtmosphereElevationChange(self.grid, model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check that the temperature and precipitation changed accordingly
        modifier.update(self.geometry, 0, 1)

        C = config.get_number("atmosphere.precip_exponential_factor_for_temperature")
        dT = -self.precip_dTdz * self.dz / 1000.0
        P = sample(model.precipitation())
        dP = np.exp(C * dT) * P - P

        check_modifier(model, modifier, T=self.dT, P=dP,
                       ts=[0.5], Ts=[self.dT], Ps=[dP])
