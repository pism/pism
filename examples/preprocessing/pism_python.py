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
# \verbatim $ pism -i foo.nc -bootstrap -Mx 41 -My 41 -Mz 21 -Lz 4000 -Mbz 5 -Lbz 500 -y 1 \endverbatim

import sys
import time
import numpy as np

try:
    import xarray as xr
except ImportError:
    print("xarray is not installed!")
    sys.exit(1)

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

# Build the dataset (dimension transpose is standard: "float thk(y, x)").
ds = xr.Dataset(
    coords={
        "x": ("x", x.astype("f4"), {
            "units": "m", "long_name": "easting",
            "standard_name": "projection_x_coordinate",
        }),
        "y": ("y", y.astype("f4"), {
            "units": "m", "long_name": "northing",
            "standard_name": "projection_y_coordinate",
        }),
    },
    data_vars={
        "topg": (("y", "x"), topg.astype("f4"),
                 {"units": "m", "standard_name": "bedrock_altitude"}),
        "thk": (("y", "x"), thk.astype("f4"),
                {"units": "m", "standard_name": "land_ice_thickness"}),
        "climatic_mass_balance": (("y", "x"), acab.astype("f4"),
                                  {"units": "m year-1",
                                   "standard_name": "land_ice_surface_specific_mass_balance"}),
        "ice_surface_temp": (("y", "x"), artm.astype("f4"),
                             {"units": "kelvin"}),
    },
    attrs={
        "Conventions": "CF-1.4",
        "history": time.asctime() + ': ' + ' '.join(sys.argv) + '\n',
    },
)

encoding = {name: {"_FillValue": np.float32(np.nan)}
            for name in ("topg", "thk", "climatic_mass_balance",
                         "ice_surface_temp")}
ds.to_netcdf(ncfile, format="NETCDF3_CLASSIC", encoding=encoding)

print('  PISM-bootable NetCDF file %s written' % ncfile)
print('  for example, run:')
print('    $ pism -i foo.nc -bootstrap -Mx 41 -My 41 -Mz 21 -Lz 4000 -Mbz 5 -Lbz 500 -y 1')
