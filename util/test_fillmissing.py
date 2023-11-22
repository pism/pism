# Copyright (C) 2023 Andy Aschwanden
#
# This file is part of pism.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import numpy as np
from numpy.testing import assert_array_almost_equal
import pandas as pd
import pytest
from typing import Tuple
import xarray as xr

from .fill_missing import laplace

@pytest.fixture(name="test_data")
def fixture_test_data() -> Tuple[xr.Dataset, xr.Dataset]:
    """
    Return a dataset with missing values
    """

    start_date = "1980-01-01"
    calendar = "standard"
    periodicity = "YS"
    reference_time = pd.Timestamp("1980-01-01")

    nx = 101
    ny = 101
    nt = 2

    dates = xr.cftime_range(start_date, periods=nt, freq=periodicity, calendar=calendar)
    x = np.linspace(-100_000, 100_000, nx)
    y = np.linspace(-100_000, 100_000, ny)
    t = np.arange(nt)
    X, Y = np.meshgrid(x, y)

    t_mask = ((np.abs(X)<25_000) & (np.abs(Y)<25_000)).reshape(1, ny, nx)
    t_mask = np.repeat(t_mask, repeats=nt, axis=0)
    temperature = np.exp(((X/100_000)**2 + (Y/100_000)**2)).reshape(1, ny, nx)
    temperature = np.repeat(temperature, repeats=nt, axis=0)

    p_mask = ((np.abs(X)<25_000) & (Y<-25_000)).reshape(1, ny, nx)
    p_mask = np.repeat(p_mask, repeats=nt, axis=0)
    precipitation = (np.exp(((X/100_000)**2 + (Y/100_000)**2)) ** 2).reshape(1, ny, nx)
    precipitation = np.repeat(precipitation, repeats=nt, axis=0)

    ds_true = xr.Dataset(
        data_vars=dict(
            temperature=(["time", "y", "x"], temperature, {"units": "Celsius"}),
            precipitation=(["time", "y", "x"], precipitation, {"units": "mm/day"}),
        ),
        coords=dict(
            time=dates,
            y=(["y"], y, {"_FillValue": False, "units": "m", "axis": "Y", "standard_name": "projection_y_coordinate"}),
            x=(["x"], x, {"_FillValue": False, "units": "m", "axis": "X", "standard_name": "projection_x_coordinate"}),
            reference_time=reference_time,
        ),
        attrs=dict(description="Test data."),
    )

    temperature = np.ma.array(data=temperature, mask=t_mask)
    precipitation = np.ma.array(precipitation, mask=p_mask)
    ds_masked = xr.Dataset(
        data_vars=dict(
            temperature=(["time", "y", "x"], temperature, {"units": "Celsius"}),
            precipitation=(["time", "y", "x"], precipitation, {"units": "mm/day"}),
        ),
        coords=dict(
            time=dates,
            y=(["y"], y, {"_FillValue": False, "units": "m", "axis": "Y", "standard_name": "projection_y_coordinate"}),
            x=(["x"], x, {"_FillValue": False, "units": "m", "axis": "X", "standard_name": "projection_x_coordinate"}),
            reference_time=reference_time,
        ),
        attrs=dict(description="Test data."),
    )
    return (ds_true, ds_masked)

def test_fill_missing(test_data):
    ds_true, ds_masked = test_data

    for t in range(len(ds_masked["time"])):
        for m_var in ["temperature"]:
            print(t, m_var)
            arr = ds_masked[m_var][t, ...].to_masked_array()
            data = arr.data
            mask = arr.mask
            
            laplace(data, mask, -1, 1.0e-2, initial_guess="mean")
            data_true = ds_true[m_var][t, ...].to_numpy()
            assert_array_almost_equal(data, data_true, decimal=1)
