#!/usr/bin/env python3

"Test the 'Given' bed deformation model."

import PISM
from PISM.testing import shallow_grid, create_forcing
import os

from unittest import TestCase

ctx = PISM.Context()
ctx.log.set_threshold(1)

class BeddefGiven(TestCase):
    def setUp(self):

        self.geometry_filename = "beddef_given_input.nc"
        self.filename          = "beddef_given.nc"
        self.ref_filename      = "beddef_given_reference.nc"

        self.grid = shallow_grid(Mx=3, My=3)

        # Create time-dependent topography changes
        create_forcing(self.grid, self.filename, "topg_delta", "meters",
                       values=[2, -4, 3], times=[1, 2, 3])

        # Create the reference topography
        topg_ref = PISM.IceModelVec2S(self.grid, "topg", PISM.WITHOUT_GHOSTS)
        topg_ref.set_attrs("bed", "bed elevation change relative to reference", "m", "m", "", 0)
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

        model.update(geometry.ice_thickness,
                     geometry.sea_level_elevation,
                     1, 1)

        topg_0 = model.bed_elevation().numpy()[0, 0]

        model.update(geometry.ice_thickness,
                     geometry.sea_level_elevation,
                     2, 1)

        topg_1 = model.bed_elevation().numpy()[0, 0]

        # -4 - 2 == -6 (see the create_forcing() call above)
        assert topg_1 - topg_0 == -6
