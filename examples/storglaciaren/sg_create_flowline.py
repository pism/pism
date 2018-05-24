#!/usr/bin/env python
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

import numpy as np

try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

x, topg, thk = np.loadtxt('sg_35m_flowline.txt', unpack=True)

output = 'storglaciaren_flowline.nc'

# Write the data:
print("Writing the data to '%s'... " % output)
nc = CDF(output, "w")
nc.createDimension("x", size=len(x))

x_var = nc.createVariable("x", 'f', dimensions=("x",))
x_var.units = "m"
x_var[:] = x

topg_var = nc.createVariable("topg", 'f', dimensions=("x",))
topg_var.units = "m"
topg_var.standard_name = "bedrock_altitude"
topg_var[:] = topg

thk_var = nc.createVariable("thk", 'f', dimensions=("x",))
thk_var.units = "m"
thk_var.standard_name = "land_ice_thickness"
thk_var[:] = thk

usurf_var = nc.createVariable("usurf", 'f', dimensions=("x",))
usurf_var.units = "m"
usurf_var.standard_name = "surface_altitude"
usurf_var[:] = topg + thk

qgeo = 0.042
bheatflx_var = nc.createVariable("bheatflx", 'f', dimensions=("x",))
bheatflx_var.units = "W m-2"
bheatflx_var[:] = qgeo * np.ones_like(x)

# generate (somewhat) reasonable acab
acab_max = 2.5  # m/a
acab_min = -3.0  # m/a
np.linspace(acab_max, acab_min, 100)
acab = np.ones_like(x)
acab[:5] = 0
acab[5:105] = np.linspace(acab_max, acab_min, 100)
acab[105:] = acab_min

ice_density = 910.0             # kg m-3
acab_var = nc.createVariable("climatic_mass_balance", 'f', dimensions=("x",))
acab_var.units = "kg m-2 year-1"
acab_var.standard_name = "land_ice_surface_specific_mass_balance_flux"
# convert from m/year ice equivalent into [kg m-2 year-1]:
acab_var[:] = acab * ice_density

# Set boundary conditions
# ------------------------------------------------------------------------------
#
# (A) Surface temperature for temperature equation bc
Tma = -6.0  # degC, mean annual air temperature at Tarfala
zcts = 1400   # m a.s.l.; altitude where CTS is at the surface


artm = np.zeros_like(x)
artm[(topg + thk) < zcts] = Tma
artm_var = nc.createVariable("ice_surface_temp", 'f', dimensions=("x",))
artm_var.units = "deg_C"
artm_var[:] = artm

nc.close()
