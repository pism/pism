#!/usr/bin/env python
"""
Tests of PISM's atmosphere models and modifiers.
"""

import PISM
from PISM.testing import *
import os
import numpy as np
from unittest import TestCase

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

options = PISM.PETSc.Options()

def check_modifier(model, modifier, T=0.0, P=0.0):
    check_difference(modifier.mean_precipitation(),
                     model.mean_precipitation(),
                     P)

    check_difference(modifier.mean_annual_temp(),
                     model.mean_annual_temp(),
                     T)

    # FIXME: add checks for point-wise time series

class PIK(TestCase):
    """Test that all the code in atmosphere::PIK runs. Does not check computed values."""
    def setUp(self):
        self.filename = "atmosphere_pik_input.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)

        precip = PISM.IceModelVec2S(self.grid, "precipitation", PISM.WITHOUT_GHOSTS)
        precip.set_attrs("climate", "dummy precipitation field", "kg m-2 s-1", "")
        precip.set(10.0)
        precip.dump(self.filename)

        config.set_string("atmosphere.pik.file", self.filename)

    def test_atmosphere_pik(self):
        "Model PIK"

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

class DeltaT(TestCase):
    def setUp(self):
        self.filename = "delta_T_input.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dT = -5.0

        create_scalar_forcing(self.filename, "delta_T", "Kelvin", [self.dT], [0])

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_delta_t(self):
        "Modifier Delta_T"

        modifier = PISM.AtmosphereDeltaT(self.grid, self.model)

        options.setValue("-atmosphere_delta_T_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier, T=self.dT, P=0.0)

class DeltaP(TestCase):
    def setUp(self):
        self.filename = "delta_P_input.nc"
        self.grid = shallow_grid()
        self.geometry = PISM.Geometry(self.grid)
        self.geometry.ice_thickness.set(1000.0)
        self.model = PISM.AtmosphereUniform(self.grid)
        self.dP = 5.0

        create_scalar_forcing(self.filename, "delta_P", "kg m-2 year-1", [self.dP], [0])

    def tearDown(self):
        os.remove(self.filename)

    def test_atmosphere_delta_p(self):
        "Modifier Delta_P"

        modifier = PISM.AtmosphereDeltaP(self.grid, self.model)

        options.setValue("-atmosphere_delta_P_file", self.filename)

        modifier.init(self.geometry)
        modifier.update(self.geometry, 0, 1)

        check_modifier(self.model, modifier,
                       P=PISM.util.convert(self.dP, "kg m-2 year-1", "kg m-2 s-1"),
                       T=0.0)

def test_atmosphere_given():
    pass

def test_atmosphere_searise_greenland():
    pass

def test_atmosphere_yearly_cycle():
    pass

def test_atmosphere_one_station():
    pass

def test_atmosphere_uniform():
    pass

def test_atmosphere_anomaly():
    pass

def test_atmosphere_paleo_precip():
    pass

def test_atmosphere_frac_p():
    pass

def test_atmosphere_lapse_rate():
    pass
