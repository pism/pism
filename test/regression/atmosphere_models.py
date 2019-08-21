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
config.set_double("grid.My", 5) # non-square grid
config.set_double("grid.Mz", 2)

seconds_per_year = 365 * 86400
# ensure that this is the correct year length
config.set_string("time.calendar", "365_day")

# silence models' initialization messages
PISM.Context().log.set_threshold(1)

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

def dummy_grid():
    "Create a dummy grid"
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)

class PIK(TestCase):
    """Test that all the code in atmosphere::PIK runs. Does not check computed values."""
    def setUp(self):
        self.filename = "atmosphere_pik_input.nc"
        self.grid = dummy_grid()
        self.geometry = create_geometry(self.grid)

        precip = PISM.IceModelVec2S(self.grid, "precipitation", PISM.WITHOUT_GHOSTS)
        precip.set_attrs("climate", "dummy precipitation field", "kg m-2 s-1", "")
        precip.set(10.0)
        precip.dump(self.filename)

        config.set_string("atmosphere.pik.file", self.filename)

    def runTest(self):
        "Atmosphere model PIK"

        parameterizations = ["martin",
                             "huybrechts_dewolde",
                             "martin_huybrechts_dewolde",
                             "era_interim",
                             "era_interim_sin",
                             "era_interim_lon"]

        for p in parameterizations:
            print("Testing parameterization {}...".format(p))

            config.set_string("atmosphere.pik.parameterization", p)

            model = PISM.AtmospherePIK(self.grid)
            model.init(self.geometry)

            # t and dt are irrelevant here
            model.update(self.geometry, 0, 1)

            model.init_timeseries([0, 0.5, 1])

            try:
                model.begin_pointwise_access()
                print("temperature time series: ", model.temp_time_series(0, 0))
                print("precipitation time series: ", model.precip_time_series(0, 0))
            finally:
                model.end_pointwise_access()

            print(model.mean_annual_temp().numpy())
            print(model.mean_precipitation().numpy())

    def tearDown(self):
        os.remove(self.filename)
