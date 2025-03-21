#!/usr/bin/env python3
"""
Tests of PISM's surface models and modifiers.
"""

import PISM
from PISM.testing import *
import sys
import os
import numpy as np
from unittest import TestCase, SkipTest

from PISM.util import convert

# set run duration to 1 second so that all forcing used here spans the duration of the run
time = PISM.Context().time
time.set_start(0)
time.set_end(1)

config = PISM.Context().config

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

log = PISM.Context().log
# silence models' initialization messages
log.set_threshold(1)

options = PISM.PETSc.Options()

def write_state(model, filename):
    "Write the state of the model to a file"

    PISM.util.prepare_output(filename)
    f = PISM.File(model.grid().ctx().com(),
                  filename,
                  PISM.PISM_NETCDF3,
                  PISM.PISM_READWRITE)
    model.define_model_state(f)
    model.write_model_state(f)

    diags = model.diagnostics()
    for k in diags.keys():
        diags[k].compute().write(f)

    f.close()

def probe_interface(model):
    """Prove the interface of a surface model to check if its public methods run successfully."""
    model.accumulation()
    model.layer_mass()
    model.layer_thickness()
    model.liquid_water_fraction()
    model.mass_flux()
    model.melt()
    model.runoff()
    model.temperature()

    model.max_timestep(0)

    model.diagnostics()

    # FIXME: this causes a memory leak
    # model.ts_diagnostics()

    model.grid()

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

def surface_simple(grid):
    return PISM.SurfaceSimple(grid, PISM.AtmosphereUniform(grid))

def climatic_mass_balance(grid, value):
    SMB = PISM.Scalar(grid, "climatic_mass_balance")
    SMB.metadata(0).long_name("surface mass balance").units("kg m^-2 s^-1").standard_name("land_ice_surface_specific_mass_balance_flux")
    SMB.set(value)
    return SMB

def ice_surface_temp(grid, value):
    temperature = PISM.Scalar(grid, "ice_surface_temp")
    temperature.metadata(0).long_name("ice temperature at the top surface").units("kelvin")
    temperature.set(value)
    return temperature

class Given(TestCase):
    def setUp(self):
        self.filename = filename("surface_given_input_")
        self.output_filename = filename("surface_given_output_")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.T = 272.15
        self.M = 1001.0

        output = PISM.util.prepare_output(self.filename)
        ice_surface_temp(self.grid, self.T).write(output)
        climatic_mass_balance(self.grid, self.M).write(output)
        output.close()

    def test_surface_given(self):
        "Model 'given'"
        atmosphere = PISM.AtmosphereUniform(self.grid)

        config.set_string("surface.given.file", self.filename)

        model = PISM.SurfaceGiven(self.grid, atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, self.T, 0.0, self.M, accumulation=self.M)

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class DeltaT(TestCase):
    def setUp(self):
        self.filename = filename("surface_delta_T_input_")
        self.output_filename = filename("surface_delta_T_output_")
        self.grid = shallow_grid()
        self.model = surface_simple(self.grid)
        self.dT = -5.0
        self.geometry = PISM.Geometry(self.grid)

        create_scalar_forcing(self.filename, "delta_T", "kelvin",
                              [self.dT], [0], time_bounds=[0, 1])

    def test_surface_delta_t(self):
        "Modifier 'delta_T'"

        config.set_string("surface.delta_T.file", self.filename)

        modifier = PISM.SurfaceDeltaT(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT)

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class ElevationChange(TestCase):
    def setUp(self):
        self.filename = filename("surface_reference_surface_")
        self.output_filename = filename("surface_lapse_rates_output_")
        self.grid     = shallow_grid()
        self.dTdz     = 1.0         # 1 kelvin per km
        self.dSMBdz   = 2.0         # m year-1 per km
        self.dz       = 1500.0      # m

        self.geometry = PISM.Geometry(self.grid)

        # save current surface elevation to use it as a "reference" surface elevation
        self.geometry.ice_surface_elevation.dump(self.filename)

        config.set_string("surface.elevation_change.file", self.filename)
        config.set_number("surface.elevation_change.temperature_lapse_rate", self.dTdz)
        config.set_number("surface.elevation_change.smb.lapse_rate", self.dSMBdz)

        self.dT = self.dz * convert(-self.dTdz, "kelvin / km", "kelvin / m")

    def test_elevation_change_shift(self):
        "Modifier elevation_change: shift"

        config.set_string("surface.elevation_change.smb.method", "shift")

        model    = surface_simple(self.grid)
        modifier = PISM.SurfaceElevationChange(self.grid, model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check changes in outputs
        modifier.update(self.geometry, 0, 1)

        ice_density = config.get_number("constants.ice.density")

        dSMB = self.dz * ice_density * convert(-self.dSMBdz,
                                               "kg m-2 year-1 / km",
                                               "kg m-2 s-1 / m")
        # shifting SMB by dSMB reduced accumulation to zero
        dA = 0.0 - sample(model.accumulation())
        dM = dA - dSMB
        dR = dM

        check_modifier(model, modifier, T=self.dT, SMB=dSMB,
                       accumulation=dA, melt=dM, runoff=dR)

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def test_elevation_change_scale(self):
        "Modifier elevation_change: scale"

        config.set_string("surface.elevation_change.smb.method", "scale")
        config.set_number("surface.elevation_change.smb.exp_factor", 1.0)

        model    = surface_simple(self.grid)
        modifier = PISM.SurfaceElevationChange(self.grid, model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check changes in outputs
        modifier.update(self.geometry, 0, 1)

        SMB = sample(model.mass_flux())
        C = config.get_number("surface.elevation_change.smb.exp_factor")

        dSMB = np.exp(C * self.dT) * SMB - SMB
        # SMB stays positive, so all SMB corresponds to accumulation
        dA = dSMB
        dM = dA - dSMB
        dR = dM

        check_modifier(model, modifier, T=self.dT, SMB=dSMB,
                       accumulation=dA, melt=dM, runoff=dR)

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class Elevation(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.output_filename = filename("surface_elevation_output_")

        # change geometry just to make this a bit more interesting
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)

        ctx = PISM.Context()

        # make a copy of the configuration database so we can re-initialize it from
        # options and then restore it
        self.config = PISM.config_from_options(ctx.com, ctx.unit_system)
        self.config.import_from(ctx.config)

    def tearDown(self):
        ctx = PISM.Context()
        ctx.config.import_from(self.config)

    def elevation_1_test(self):
        "Model 'elevation', test 1"
        model = PISM.SurfaceElevation(self.grid, PISM.AtmosphereUniform(self.grid))

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        T            = 268.15
        omega        = 0.0
        SMB          = -8.651032746943449e-05
        accumulation = 0.0
        melt         = -SMB
        runoff       = melt

        check_model(model, T, omega, SMB,
                    accumulation=accumulation, melt=melt, runoff=runoff)

        probe_interface(model)

    def elevation_2_test(self):
        "Model 'elevation', test 2"
        T_min = -5.0
        T_max = 0.0
        z_min = 1000.0
        z_ela = 1100.0
        z_max = 1500.0
        M_min = -1.0
        M_max = 5.0

        self.geometry.ice_thickness.set(0.5 * (z_min + z_max))
        self.geometry.ensure_consistency(0.0)

        options.setValue("-ice_surface_temp", "{},{},{},{}".format(T_min, T_max, z_min, z_max))
        options.setValue("-climatic_mass_balance",
                         "{},{},{},{},{}".format(M_min, M_max, z_min, z_ela, z_max))

        ctx = PISM.Context()
        PISM.set_config_from_options(ctx.unit_system, ctx.config)

        T = PISM.util.convert(0.5 * (T_min + T_max), "degree_Celsius", "kelvin")
        SMB = PISM.util.convert(1.87504, "m/year", "m/s") * config.get_number("constants.ice.density")

        model = PISM.SurfaceElevation(self.grid, PISM.AtmosphereUniform(self.grid))

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, T=T, SMB=SMB, omega=0, mass=0, thickness=0, accumulation=SMB)

class TemperatureIndex1(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.atmosphere = PISM.AtmosphereUniform(self.grid)
        self.output_filename = filename("surface_pdd_output_")

    def pdd_test(self):
        "Model 'pdd', test 1"
        config.set_string("surface.pdd.method", "expectation_integral")

        model = PISM.SurfaceTemperatureIndex(self.grid, self.atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        T = config.get_number("atmosphere.uniform.temperature", "kelvin")
        omega = 0.0

        accumulation = config.get_number("atmosphere.uniform.precipitation", "kg m-2 second-1")
        melt         = accumulation
        runoff       = melt * (1.0 - config.get_number("surface.pdd.refreeze"))
        SMB          = accumulation - runoff

        check_model(model, T, omega, SMB, accumulation=accumulation, melt=melt,
                    runoff=runoff)

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.output_filename)

class TemperatureIndex2(TestCase):
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
        config.set_number("surface.pdd.std_dev.value", 0.0)
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

        check_model(model, T=self.T, SMB=self.SMB, omega=0.0, mass=0.0, thickness=0.0,
                    melt=40, runoff=16)

class PIK(TestCase):
    def setUp(self):
        self.filename = filename("surface_pik_input_")
        self.output_filename = filename("surface_pik_output_")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.M = 1001.0
        self.T = 233.13

        self.geometry.latitude.set(-80.0)
        self.geometry.ice_thickness.set(2000.0)
        self.geometry.ensure_consistency(0.0)

        climatic_mass_balance(self.grid, self.M).dump(self.filename)

        config.set_string("input.file", self.filename)

    def surface_pik_test(self):
        "Model 'pik'"
        model = PISM.SurfacePIK(self.grid, PISM.AtmosphereUniform(self.grid))

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        check_model(model, self.T, 0.0, self.M, accumulation=self.M)

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)
        config.set_string("input.file", "")

class Simple(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.output_filename = filename("surface_simple_output_")
        self.atmosphere = PISM.AtmosphereUniform(self.grid)
        self.geometry = PISM.Geometry(self.grid)

    def simple_test(self):
        "Model 'simple'"
        atmosphere = self.atmosphere

        model = PISM.SurfaceSimple(self.grid, atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        T = sample(atmosphere.air_temperature())
        M = sample(atmosphere.precipitation())

        check_model(model, T, 0.0, M, accumulation=M)

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.output_filename)

class Anomaly(TestCase):
    def setUp(self):
        self.filename = filename("surface_anomaly_input_")
        self.output_filename = filename("surface_anomaly_output_")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.model = surface_simple(self.grid)
        self.dSMB = -(config.get_number("atmosphere.uniform.precipitation", "kg m-2 s-1") + 5.0)
        self.dT = 2.0

        PISM.util.prepare_output(self.filename)

        delta_SMB = PISM.Scalar(self.grid, "climatic_mass_balance_anomaly")
        delta_SMB.metadata(0).long_name("2D surface mass flux anomaly").units("kg m^-2 s^-1").output_units("kg m^-2 s^-1")
        delta_SMB.set(self.dSMB)

        delta_SMB.write(self.filename)

        delta_T = PISM.Scalar(self.grid, "ice_surface_temp_anomaly")
        delta_T.metadata(0).long_name("2D surface temperature anomaly").units("kelvin").output_units("kelvin")
        delta_T.set(self.dT)

        delta_T.write(self.filename)

    def anomaly_test(self):
        "Modifier 'anomaly'"

        config.set_string("surface.anomaly.file", self.filename)

        modifier = PISM.SurfaceAnomaly(self.grid, self.model)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        # once anomaly is applied the SMB is negative, so the new accumulation is zero
        dA = 0.0 - sample(self.model.accumulation())
        dM = dA - self.dSMB
        dR = dM

        check_modifier(self.model, modifier, T=self.dT, SMB=self.dSMB,
                       accumulation=dA, melt=dM, runoff=dR)

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class Cache(TestCase):
    def setUp(self):
        self.filename = filename("surface_dT_")
        self.output_filename = filename("surface_cache_output_")
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.simple = surface_simple(self.grid)

        time_bounds = np.array([0, 1, 1, 2, 2, 3, 3, 4]) * seconds_per_year
        create_scalar_forcing(self.filename, "delta_T", "kelvin", [1, 2, 3, 4],
                              times=None, time_bounds=time_bounds)

        config.set_string("surface.delta_T.file", self.filename)

        self.delta_T = PISM.SurfaceDeltaT(self.grid, self.simple)

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

            # the "cache" modifier will (note: times in years)
            # - timestep [0, 1]: call update(0, 1), which evaluates forcing at 0.5.
            #
            #   Since we assume that times are at midpoints of intervals indicated using
            #   bounds, this will return 1.
            #
            # - timestep [1, 2]: re-use cached values, i.e. return 1
            #
            # - timestep [2, 3]: call update(2, 1), which evaluates forcing at 2.5
            #
            #   Since we assume that times are at midpoints of intervals indicated using
            #   bounds, this will return 3.
            #
            # - timestep [3, 4]: re-use cached values, i.e. return 3.

            diff.append(cached - original)

        write_state(modifier, self.output_filename)

        probe_interface(modifier)

        np.testing.assert_almost_equal(diff, [1, 1, 3, 3])

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class ForceThickness(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.model = surface_simple(self.grid)
        self.filename = filename("surface_force_to_thickness_input_")
        self.output_filename = filename("surface_force_to_thickness_output_")

        self.H = 1000.0
        self.dH = 1000.0

        self.geometry.ice_thickness.set(self.H)

        # save ice thickness to a file to use as the target thickness
        PISM.util.prepare_output(self.filename)
        self.geometry.ice_thickness.write(self.filename)

        ftt_mask = PISM.Scalar(self.grid, "ftt_mask")
        ftt_mask.set(1.0)
        ftt_mask.write(self.filename)

        alpha       = 10.0
        ice_density = config.get_number("constants.ice.density")
        self.dSMB   = -ice_density * alpha * self.dH

        config.set_string("surface.force_to_thickness.file", self.filename)
        config.set_number("surface.force_to_thickness.alpha", convert(alpha, "1/s", "1/year"))

    def forcing_test(self):
        "Modifier ForceThickness"
        modifier = PISM.SurfaceForceThickness(self.grid, self.model)

        modifier.init(self.geometry)

        self.geometry.ice_thickness.set(self.H + self.dH)

        modifier.update(self.geometry, 0, 1)

        dA   = 0.0 - sample(self.model.accumulation())
        dM   = dA - self.dSMB
        dR   = dM

        check_modifier(self.model, modifier, SMB=self.dSMB,
                       accumulation=dA, melt=dM, runoff=dR)

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class EISMINTII(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.output_filename = filename("surface_eismint_output_")

    def eismintii_test(self):
        "Model EISMINTII: define and write model state; get diagnostics"

        for experiment in "ABCDEFGHIJKL":
            model = PISM.SurfaceEISMINTII(self.grid, ord(experiment))

            model.init(self.geometry)

            model.update(self.geometry, 0, 1)

            write_state(model, self.output_filename)
            probe_interface(model)

            os.remove(self.output_filename)

    def tearDown(self):
        try:
            os.remove(self.output_filename)
        except:
            pass

class Initialization(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.output_filename = filename("surface_init_output_")
        self.model = surface_simple(self.grid)

    def initialization_test(self):
        "Modifier InitializationHelper"

        modifier = PISM.SurfaceInitialization(self.grid, self.model)

        modifier.init(self.geometry)

        modifier.update(self.geometry, 0, 1)

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.output_filename)

class Factory(TestCase):
    def setUp(self):
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

    def factory_test(self):
        "Surface model factory"
        atmosphere = PISM.AtmosphereUniform(self.grid)

        factory = PISM.SurfaceFactory(self.grid, atmosphere)

        simple = factory.create("simple")

        model = factory.create("simple,cache")

        try:
            factory.create("invalid_model")
            return False
        except RuntimeError:
            pass

        try:
            factory.create("simple,invalid_modifier")
            return False
        except RuntimeError:
            pass

    def tearDown(self):
        pass


class ISMIP6(TestCase):
    def prepare_reference_data(self, grid, filename):

        usurf   = PISM.model.createIceSurfaceVec(grid)
        SMB_ref = PISM.Scalar(grid, "climatic_mass_balance")
        T_ref   = PISM.Scalar(grid, "ice_surface_temp")

        usurf.metadata(0).set_string("units", "m")

        SMB_ref.metadata(0).long_name("reference SMB").units("kg m^-2 s^-1").output_units("kg m^-2 s^-1").standard_name("land_ice_surface_specific_mass_balance_flux")

        T_ref.metadata(0).set_string("units", "kelvin")

        out = PISM.util.prepare_output(filename, append_time=True)

        T_ref.set(260.0)
        SMB_ref.set(0.0)
        usurf.set(0.0)

        # write time-independent fields
        for v in [usurf, SMB_ref, T_ref]:
            v.metadata().set_time_independent(True)
            v.write(out)

        out.close()

    def prepare_climate_forcing(self, grid, filename):

        aSMB = PISM.Scalar(grid, "climatic_mass_balance_anomaly")
        aSMB.metadata(0).long_name("SMB anomaly").units("kg m^-2 s^-1").output_units("kg m^-2 s^-1")

        dSMBdz = PISM.Scalar(grid, "climatic_mass_balance_gradient")
        dSMBdz.metadata(0).long_name("SMB gradient").units("kg m^-2 s^-1 m^-1").output_units("kg m^-2 s^-1 m^-1")

        aT = PISM.Scalar(grid, "ice_surface_temp_anomaly")
        aT.metadata(0).long_name("temperature anomaly").units("kelvin").output_units("kelvin")

        dTdz = PISM.Scalar(grid, "ice_surface_temp_gradient")
        dTdz.metadata(0).long_name("surface temperature gradient").units("K m^-1").output_units("K m^-1")

        out = PISM.util.prepare_output(filename, append_time=False)

        bounds = PISM.VariableMetadata("time_bounds", self.ctx.unit_system)

        PISM.define_time_bounds(bounds, "time", "nv", out, PISM.PISM_DOUBLE)

        SMB_anomaly  = 1.0
        T_anomaly    = 1.0
        SMB_gradient = 1.0
        T_gradient   = 1.0

        # monthly steps
        dt = (365 * 86400) / 12.0

        for j in range(12):

            t = self.ctx.time.current() + j * dt

            PISM.append_time(out, self.ctx.config, t)

            PISM.write_time_bounds(out, bounds, j, [t, t + dt])

            aSMB.set(t * SMB_anomaly)

            dSMBdz.set(t * SMB_gradient)

            aT.set(t * T_anomaly)

            dTdz.set(t * T_gradient)

            for v in [aSMB, dSMBdz, aT, dTdz]:
                v.write(out)

        out.redef()
        out.write_attribute("time", "bounds", "time_bounds")
        out.write_attribute("time", "units", "seconds since 1-1-1")

        out.close()

    def setUp(self):

        self.ctx = PISM.Context()

        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        self.geometry.ice_surface_elevation.set(100.0)

        self.forcing_file = filename("surface_ismip6_forcing_")
        self.reference_file = filename("surface_ismip6_reference_")

        self.prepare_reference_data(self.grid, self.reference_file)
        self.prepare_climate_forcing(self.grid, self.forcing_file)

        self.ctx.config.set_string("surface.ismip6.file", self.forcing_file)
        self.ctx.config.set_string("surface.ismip6.reference_file", self.reference_file)

    def ismip6_test(self):
        "Surface model ISMIP6"

        atmosphere = PISM.AtmosphereUniform(self.grid)

        model = PISM.SurfaceISMIP6(self.grid, atmosphere)

        model.init(self.geometry)

        t = self.ctx.time.current()
        dt = model.max_timestep(t).value()

        model.update(self.geometry, t, dt)

    def tearDown(self):
        os.remove(self.reference_file)
        os.remove(self.forcing_file)

if __name__ == "__main__":

    t = ISMIP6()

    t.setUp()
    t.ismip6_test()
    t.tearDown()
