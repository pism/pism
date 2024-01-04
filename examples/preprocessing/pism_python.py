#!/usr/bin/env python3

# Copyright (C) 2009-2015, 2018, 2024 the PISM Authors

# @package pism_python
# \author the PISM authors
# \brief Creates "from scratch" a boring dataset with the right format
# to use as a PISM bootstrapping file.
# \details Example use of Python for this purpose.
#
# Usage, including a minimal PISM call to bootstrap from this file:
#
# \verbatim $ pism_python.py  # creates foo.nc \endverbatim
# \verbatim $ pismr -i foo.nc -bootstrap -Mx 41 -My 41 -Mz 21 -Lz 4000 -Mbz 5 -Lbz 500 -y 1 \endverbatim

import sys
import time
import numpy as np

import xarray as xr

# set up the grid:
Lx = 1e6
Ly = 1e6
Mx = 51
My = 71
x = np.linspace(-Lx, Lx, Mx)
y = np.linspace(-Ly, Ly, My)

# create dummy fields
[xx, yy] = np.meshgrid(x, y)  # if there were "ndgrid" in numpy we would use it
acab = np.zeros((My, Mx))
artm = np.zeros((My, Mx)) + 273.15 + 10.0  # 10 degrees Celsius
topg = 1000.0 + 200.0 * (xx + yy) / max(Lx, Ly)  # change "1000.0" to "0.0" to test
# flotation criterion, etc.
thk = 3000.0 * (1.0 - 3.0 * (xx ** 2 + yy ** 2) / Lx ** 2)
thk[thk < 0.0] = 0.0

# Output filename
ncfile = 'foo.nc'

historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'

coords = {
    "x": (
        ["x"],
        x,
        {
            "units": "m",
            "axis": "X",
            "standard_name": "projection_x_coordinate",
            "long_name": "x-coordinate in projected coordinate system",
        },
    ),
    "y": (
        ["y"],
        y,
        {
            "units": "m",
            "axis": "Y",
            "standard_name": "projection_y_coordinate",
            "long_name": "y-coordinate in projected coordinate system",
        },
    ),
}

ds = xr.Dataset(
    {
        "topg": xr.DataArray(
            data=topg,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "standard_name": "bedrock_altitude",
                "units": "m",
            },
        ),
        "thk": xr.DataArray(
            data=thk,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "standard_name": "land_ice_thickness",
                "units": "m",
            },
        ),
        "ice_surface_temp": xr.DataArray(
            data=artm,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "units": "K",
            },
        ),
        "climatic_mass_balance": xr.DataArray(
            data=acab,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "standard_name": "land_ice_surface_specific_mass_balance",
                "units": "m year-1",
            },
        ),
    },
    attrs={"Conventions": "CF-1.7", "history": historystr},
)
ds.to_netcdf(ncfile)
