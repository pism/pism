#!/usr/bin/env python3
"""
Tests of PISM's debris models and of the debris-related ice melt enhancement.
"""

import PISM
import os
import numpy
from unittest import TestCase
from PISM.testing import filename as tmp_name

config = PISM.Context().config

# reduce the grid size to speed this up
config.set_number("grid.Mx", 3)
config.set_number("grid.My", 3)
config.set_number("grid.Mz", 5)
config.set_number("grid.Lx", 1e2)  # in km
config.set_number("grid.Ly", 1e2)  # in km

log = PISM.Context().log
# silence models' initialization messages
log.set_threshold(1)

# debris thicknesses (meters) covering all three branches of equation (14)
THICKNESS = [0.0, 0.005, 0.015, 0.03, 0.041, 0.1, 0.5, 1.5, 3.0]

# corresponding values of the melt factor
FACTOR = [1.0,
          26.667 * 0.005 + 1.0,
          26.667 * 0.015 + 1.0,
          -16.0 * 0.03 + 1.64,
          -16.0 * 0.041 + 1.64,
          0.1061 * 0.1 ** -0.7205 - 0.07922,
          0.1061 * 0.5 ** -0.7205 - 0.07922,
          1e-3,
          1e-3]


def create_grid():
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.config, ctx.size)
    return PISM.Grid(ctx.ctx, params)


def create_geometry(grid):
    geometry = PISM.Geometry(grid)

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ice_thickness.set(1000.0)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.ensure_consistency(0.0)

    return geometry


def create_input_file(filename, grid, thickness, melt_factor):
    "Write 'debris_thickness' and 'debris_melt_factor' to a file."
    output = PISM.util.prepare_output(filename)

    H = PISM.Scalar(grid, "debris_thickness")
    H.metadata(0).long_name("debris thickness").units("m")
    H.set(thickness)
    output.define_variable(H.metadata())
    H.write(output)

    F = PISM.Scalar(grid, "debris_melt_factor")
    F.metadata(0).long_name("sub-debris melt enhancement factor").units("1")
    F.set(melt_factor)
    output.define_variable(F.metadata())
    F.write(output)

    output.close()


def sample(vec):
    return vec.to_numpy()[0, 0]


class Equation14(TestCase):
    def test_verhaegen_melt_factor(self):
        "Equation (14) of Verhaegen and Huybrechts (2026)"

        f = PISM.IceMeltEnhancement.verhaegen_melt_factor

        for h, expected in zip(THICKNESS, FACTOR):
            numpy.testing.assert_almost_equal(f(h), expected)

        # the factor peaks at the effective debris thickness (1.5 cm)
        numpy.testing.assert_almost_equal(f(0.015), 1.4, decimal=4)

        # melt is enhanced for thin debris and suppressed for thick debris
        self.assertGreater(f(0.005), 1.0)
        self.assertLess(f(0.1), 1.0)

        # the factor never drops below 1e-3
        self.assertEqual(f(1e6), 1e-3)

        # negative thicknesses (and NaNs) are treated as "no debris"
        self.assertEqual(f(-1.0), 1.0)
        self.assertEqual(f(float("nan")), 1.0)


class IceMeltEnhancement(TestCase):
    def setUp(self):
        self.filename = tmp_name("debris_input")
        self.grid = create_grid()
        self.geometry = create_geometry(self.grid)

        self.thickness = 0.1
        self.melt_factor = 0.25

        create_input_file(self.filename, self.grid, self.thickness, self.melt_factor)

        config.set_string("debris.given.file", self.filename)
        config.set_string("debris.ice_melt_enhancement.file", self.filename)

    def tearDown(self):
        os.remove(self.filename)
        config.set_string("debris.ice_melt_enhancement.model", "none")
        config.set_string("debris.given.file", "")
        config.set_string("debris.ice_melt_enhancement.file", "")

    def test_none(self):
        "Model 'none': the factor is 1 everywhere"

        config.set_string("debris.ice_melt_enhancement.model", "none")

        model = PISM.IceMeltEnhancement(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        numpy.testing.assert_almost_equal(sample(model.ice_melt_enhancement()), 1.0)
        self.assertTrue(model.max_timestep(0, None).infinite())

    def test_given(self):
        "Model 'given': read 'debris_melt_factor' from a file"

        config.set_string("debris.ice_melt_enhancement.model", "given")

        model = PISM.IceMeltEnhancement(self.grid)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        numpy.testing.assert_almost_equal(sample(model.ice_melt_enhancement()),
                                          self.melt_factor)

    def test_verhaegen(self):
        "Model 'verhaegen': apply equation (14) to the debris thickness"

        config.set_string("debris.ice_melt_enhancement.model", "verhaegen")

        debris = PISM.DebrisGiven(self.grid)
        model = PISM.IceMeltEnhancement(self.grid, debris)
        model.init(self.geometry)
        model.update(self.geometry, 0, 1)

        numpy.testing.assert_almost_equal(sample(debris.debris()), self.thickness)
        numpy.testing.assert_almost_equal(
            sample(model.ice_melt_enhancement()),
            PISM.IceMeltEnhancement.verhaegen_melt_factor(self.thickness))

    def test_verhaegen_requires_a_debris_model(self):
        "Model 'verhaegen' fails without a debris model"

        config.set_string("debris.ice_melt_enhancement.model", "verhaegen")

        with self.assertRaises(RuntimeError):
            PISM.IceMeltEnhancement(self.grid)

    def test_invalid_model(self):
        "An unsupported model name is an error"

        config.set_string("debris.ice_melt_enhancement.model", "invalid")

        with self.assertRaises(RuntimeError):
            PISM.IceMeltEnhancement(self.grid)
