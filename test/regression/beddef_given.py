#!/usr/bin/env python3

"Test the 'Given' bed deformation model."

import PISM
from PISM.testing import shallow_grid, create_forcing
import os
import numpy as np

from unittest import TestCase

ctx = PISM.Context()
ctx.log.set_threshold(1)

times = [0.25, 0.75, 1.25, 1.75, 2.25, 2.75]

# set run duration so that all forcing used here spans the duration of the run
time = ctx.time
time.set_start(times[0])
time.set(times[0])
time.set_end(times[-1])

class BeddefGiven(TestCase):
    def setUp(self):

        self.geometry_filename = "beddef_given_input.nc"
        self.filename          = "beddef_given.nc"
        self.ref_filename      = "beddef_given_reference.nc"

        self.grid = shallow_grid(Mx=3, My=3)

        # Create time-dependent topography changes
        #
        # Note that the array `values` below contains intervals where topg_delta is
        # constant in time. We need this to get easily predictable outputs: PISM uses
        # piecewise-linear interpolation in time to compute averages of topg_delta over
        # time steps.
        create_forcing(self.grid, self.filename, "topg_delta", "meters",
                       values=[2, 2, -4, -4, 3, 3], times=times,
                       time_bounds=[0, 0.5, 1, 1.5, 2, 2.5, 3])

        # Create the reference topography
        topg_ref = PISM.Scalar(self.grid, "topg")
        topg_ref.metadata(0).long_name("reference bed elevation").units("m")
        topg_ref.set(-1.0)

        topg_ref.dump(self.ref_filename)

        # Create the "geometry" file to bootstrap from. Values set below do not matter.
        geometry = PISM.Geometry(self.grid)
        geometry.ice_thickness.set(0.0)
        geometry.sea_level_elevation.set(0.0)
        geometry.bed_elevation.set(0.0)
        geometry.ensure_consistency(0.0)
        geometry.dump(self.geometry_filename)

        # Set configuration flags
        config = self.grid.ctx().config()

        config.set_flag("input.bootstrap", True)
        config.set_string("input.file", self.geometry_filename)
        config.set_string("bed_deformation.given.file", self.filename)
        config.set_string("bed_deformation.given.reference_file", self.ref_filename)
        config.set_number("bed_deformation.update_interval", 0.0)

    def tearDown(self):
        try:
            # remove files
            os.remove(self.geometry_filename)
            os.remove(self.filename)
            os.remove(self.ref_filename)
        except:
            pass

    def bed_def_given_test(self):
        "Test -bed_def given"

        model = PISM.GivenTopography(self.grid)

        geometry = PISM.Geometry(self.grid)
        geometry.ice_thickness.set(0.0)
        geometry.sea_level_elevation.set(0.0)

        opts = PISM.process_input_options(ctx.com, ctx.config)
        model.init(opts, geometry.ice_thickness, geometry.sea_level_elevation)

        # step from t=0.25 to t=0.75
        model.update(geometry.ice_thickness,
                     geometry.sea_level_elevation,
                     0.25, 0.5)

        topg_0 = model.bed_elevation().to_numpy()[0, 0]

        # step from t=0.75 to t=1.25
        #
        # this update is used to advance the time used by the bed deformation model
        model.update(geometry.ice_thickness,
                     geometry.sea_level_elevation,
                     0.75, 0.5)

        # step from t=1.25 to t=1.75
        model.update(geometry.ice_thickness,
                     geometry.sea_level_elevation,
                     1.25, 0.5)

        topg_1 = model.bed_elevation().to_numpy()[0, 0]

        # -4 - 2 == -6 (see the create_forcing() call above)
        np.testing.assert_almost_equal(topg_1 - topg_0, -6)
