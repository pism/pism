#!/usr/bin/env python
"""
Tests of PISM's atmosphere models and modifiers.
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
config.set_double("grid.My", 6) # non-square grid
config.set_double("grid.Mz", 5)

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

# silence models' initialization messages
log = PISM.Context().log.set_threshold(1)

def create_geometry(grid):
    geometry = PISM.Geometry(grid)

    geometry.cell_area.set(grid.dx() * grid.dy())
    geometry.latitude.set(0.0)
    geometry.longitude.set(0.0)

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)

    geometry.ice_thickness.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)

    geometry.ensure_consistency(0.0)

    return geometry

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
    "Check if vec[0,0] is almost equal to value."
    np.testing.assert_almost_equal(vec.numpy()[0,0], value)

def check_difference(A, B, value):
    "Check if the difference between A and B is almost equal to value."
    delta = A.numpy() - B.numpy()
    np.testing.assert_almost_equal(delta[0,0], value)

def check_model(model, T, SMB):
    check(model.mean_precipitation(), SMB)
    check(model.mean_annual_temp(), T)

def check_modifier(model, modifier, T=0.0, P=0.0):
    check_difference(modifier.mean_precipitation(),
                     model.mean_precipitation(),
                     P)

    check_difference(modifier.mean_annual_temp(),
                     model.mean_annual_temp(),
                     T)

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
        self.grid = dummy_grid()
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dT = -5.0
        self.geometry = create_geometry(self.grid)

        create_dummy_forcing_file(self.filename, "delta_T", "Kelvin", self.dT)

    def runTest(self):
        "Modifier Delta_T"


        modifier = PISM.AtmosphereDeltaT(self.grid, self.model)

        options.setValue("-surface_delta_T_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT)

    def tearDown(self):
        os.remove(self.filename)

class ElevationCorrection(TestCase):
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

        config.set_string("surface.elevation_correction.file", self.filename)
        config.set_double("surface.elevation_correction.temperature.lapse_rate", self.dTdz)

    def runTest(self):
        "Modifier lapse_rate"

        modifier = PISM.SurfaceElevationCorrection(self.grid, self.model)

        modifier.init(self.geometry)

        # change surface elevation
        self.geometry.ice_surface_elevation.shift(self.dz)

        # check that the temperature changed accordingly
        modifier.update(self.geometry, 0, 1)
        check_modifier(self.model, modifier, T=self.dT)

    def tearDown(self):
        os.remove(self.filename)


if __name__ == "__main__":

    t = ElevationCorrection()

    t.setUp()
    t.runTest()
    t.tearDown()
