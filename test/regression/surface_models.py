#!/usr/bin/env python
"""
Tests of PISM's surface models and modifiers.
"""

import PISM
from PISM.testing import *
import os
import numpy as np
from unittest import TestCase

config = PISM.Context().config

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

# silence models' initialization messages
PISM.Context().log.set_threshold(1)

options = PISM.PETSc.Options()

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

class DeltaT(TestCase):
    def setUp(self):
        self.filename = "delta_T_input.nc"
        self.grid = shallow_grid()
        self.model = surface_simple(self.grid)
        self.dT = -5.0
        self.geometry = PISM.Geometry(self.grid)

        create_scalar_forcing(self.filename, "delta_T", "Kelvin", [self.dT], [0])

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
        self.grid = shallow_grid()
        self.model = surface_simple(self.grid)
        self.dTdz = 1.0         # 1 Kelvin per km
        self.dT = -1.0
        self.dz = 1000.0

        self.geometry = PISM.Geometry(self.grid)

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

def test_surface_elevation():
    pass

def test_surface_given():
    pass

def test_surface_pdd():
    pass

def test_surface_pik():
    pass

def test_surface_simple():
    pass

def test_surface_anomaly():
    pass

def test_surface_cache():
    pass

def test_surface_forcing():
    pass
