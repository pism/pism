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

from PISM.util import convert

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

def write_state(model, filename):
    "Write the state of the model to a file"

    PISM.util.prepare_output(filename)
    f = PISM.PIO(model.grid().ctx().com(), "netcdf3",
                 filename, PISM.PISM_READWRITE)
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
        self.output_filename = "given_output.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.T = 272.15
        self.M = 1001.0

        self.create_given_input_file(self.filename, self.grid, self.T, self.M)

    def runTest(self):
        "Model given"
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
        self.filename = "delta_T_input.nc"
        self.output_filename = "delta_T_output.nc"
        self.grid = dummy_grid()
        self.model = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))
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

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class LapseRates(TestCase):
    def setUp(self):
        self.filename = "reference_surface.nc"
        self.output_filename = "lapse_rates_output.nc"
        self.grid     = dummy_grid()
        self.model    = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))
        self.dTdz     = 1.0         # 1 Kelvin per km
        self.dSMBdz   = 2.0         # m year-1 per km
        self.dz       = 1500.0      # m

        self.geometry = create_geometry(self.grid)

        # save current surface elevation to use it as a "reference" surface elevation
        self.geometry.ice_surface_elevation.dump(self.filename)

        config.set_string("surface.lapse_rate.file", self.filename)

        options.setValue("-temp_lapse_rate", self.dTdz)
        options.setValue("-smb_lapse_rate", self.dSMBdz)

        ice_density = config.get_double("constants.ice.density")

        self.dSMB = self.dz * ice_density * convert(-self.dSMBdz,
                                                    "kg m-2 year-1 / km",
                                                    "kg m-2 s-1 / m")
        self.dT = self.dz * convert(-self.dTdz, "Kelvin / km", "Kelvin / m")

    def runTest(self):
        "Modifier lapse_rate"

        modifier = PISM.SurfaceLapseRates(self.grid, self.model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check that the temperature changed accordingly
        modifier.update(self.geometry, 0, 1)

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

class Elevation(TestCase):
    def setUp(self):
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.output_filename = "elevation_output.nc"

        # change geometry just to make this a bit more interesting
        self.geometry.ice_thickness.set(1000.0)
        self.geometry.ensure_consistency(0.0)

    def runTest(self):
        "Model Elevation"
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

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.output_filename)

class TemperatureIndex(TestCase):
    def setUp(self):
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.atmosphere = PISM.AtmosphereUniform(self.grid)
        self.output_filename = "pdd_output.nc"

    def runTest(self):
        "Model TemperatureIndex"
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

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.output_filename)

class PIK(TestCase):
    def create_input(self, filename, grid, mass_flux):
        PISM.util.prepare_output(filename)

        M = PISM.IceModelVec2S(grid, "climatic_mass_balance", PISM.WITHOUT_GHOSTS)
        M.set_attrs("climate", "top surface mass balance", "kg m-2 s-1", "")
        M.set(mass_flux)
        M.write(filename)

    def setUp(self):
        self.filename = "pik_input.nc"
        self.output_filename = "pik_output.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        self.M = 1001.0
        self.T = 273.15 + 30.0

        self.create_input(self.filename, self.grid, self.M)

        config.set_string("input.file", self.filename)

    def runTest(self):
        "Model PIK"
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
        self.grid = dummy_grid()
        self.output_filename = "simple_output.nc"
        self.atmosphere = PISM.AtmosphereUniform(self.grid)
        self.geometry = create_geometry(self.grid)

    def runTest(self):
        "Model Simple"
        atmosphere = self.atmosphere

        model = PISM.SurfaceSimple(self.grid, atmosphere)

        model.init(self.geometry)

        model.update(self.geometry, 0, 1)

        T = sample(atmosphere.mean_annual_temp())
        M = sample(atmosphere.mean_precipitation())

        check_model(model, T, 0.0, M, accumulation=M)

        write_state(model, self.output_filename)
        probe_interface(model)

    def tearDown(self):
        os.remove(self.output_filename)

class Anomaly(TestCase):
    def setUp(self):
        self.filename = "surface_anomaly_input.nc"
        self.output_filename = "anomaly_output.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))
        self.dSMB = -(config.get_double("atmosphere.uniform.precipitation", "kg m-2 s-1") + 5.0)
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
        self.output_filename = "cache_output.nc"
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

        write_state(modifier, self.output_filename)
        probe_interface(modifier)

    def tearDown(self):
        os.remove(self.filename)
        os.remove(self.output_filename)

class ForceThickness(TestCase):
    def setUp(self):
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.model = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))
        self.filename = "force_to_thickness_input.nc"
        self.output_filename = "force_to_thickness_output.nc"

        self.H = 1000.0
        self.dH = 1000.0

        self.geometry.ice_thickness.set(self.H)

        # save ice thickness to a file to use as the target thickness
        PISM.util.prepare_output(self.filename)
        self.geometry.ice_thickness.write(self.filename)

        ftt_mask = PISM.IceModelVec2S(self.grid, "ftt_mask", PISM.WITHOUT_GHOSTS)
        ftt_mask.set(1.0)
        ftt_mask.write(self.filename)

        alpha       = 10.0
        ice_density = config.get_double("constants.ice.density")
        self.dSMB   = -ice_density * alpha * self.dH

        config.set_string("surface.force_to_thickness_file", self.filename)
        config.set_double("surface.force_to_thickness.alpha", convert(alpha, "1/s", "1/year"))

    def runTest(self):
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
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.output_filename = "eismint_output.nc"

    def runTest(self):
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
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)
        self.output_filename = "init_output.nc"
        self.model = PISM.SurfaceSimple(self.grid, PISM.AtmosphereUniform(self.grid))

    def runTest(self):
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
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

    def runTest(self):
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

if __name__ == "__main__":

    t = ForceThickness()

    t.setUp()
    t.runTest()
    t.tearDown()
