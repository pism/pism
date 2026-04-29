#!/usr/bin/env python3
#
# Copyright (C) 2011, 2014, 2018 Andy Aschwanden
#
# This file is part of PISM.
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

import sys
import numpy as np

try:
    import xarray as xr
except ImportError:
    print("xarray is not installed!")
    sys.exit(1)

x, topg, thk = np.loadtxt('sg_35m_flowline.txt', unpack=True)

output = 'storglaciaren_flowline.nc'

# generate (somewhat) reasonable acab
acab_max = 2.5  # m/a
acab_min = -3.0  # m/a
acab = np.ones_like(x)
acab[:5] = 0
acab[5:105] = np.linspace(acab_max, acab_min, 100)
acab[105:] = acab_min

ice_density = 910.0             # kg m-3

# Set boundary conditions
# (A) Surface temperature for temperature equation bc
Tma = -6.0  # degC, mean annual air temperature at Tarfala
zcts = 1400  # m a.s.l.; altitude where CTS is at the surface
artm = np.zeros_like(x)
artm[(topg + thk) < zcts] = Tma

qgeo = 0.042

print("Writing the data to '%s'... " % output)
ds = xr.Dataset(
    coords={"x": ("x", x.astype("f4"), {"units": "m"})},
    data_vars={
        "topg": (("x",), topg.astype("f4"),
                 {"units": "m", "standard_name": "bedrock_altitude"}),
        "thk": (("x",), thk.astype("f4"),
                {"units": "m", "standard_name": "land_ice_thickness"}),
        "usurf": (("x",), (topg + thk).astype("f4"),
                  {"units": "m", "standard_name": "surface_altitude"}),
        "bheatflx": (("x",), (qgeo * np.ones_like(x)).astype("f4"),
                     {"units": "W m^-2"}),
        "climatic_mass_balance": (("x",), (acab * ice_density).astype("f4"),
                                  {"units": "kg m^-2 year^-1",
                                   "standard_name": "land_ice_surface_specific_mass_balance_flux"}),
        "ice_surface_temp": (("x",), artm.astype("f4"), {"units": "deg_C"}),
    },
)
ds.to_netcdf(output, mode="w")
