#!/usr/bin/env python
#
# Copyright (C) 2011 Andy Aschwanden
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
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
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

x,topg,thk = np.loadtxt('sg_35m_flowline.txt',unpack=True)

output = 'storglaciaren_flowline.nc'

# Write the data:
print("Writing the data to '%s'... " % output)
nc = CDF(output, "w")
nc.createDimension("x", size=len(x))

x_var = nc.createVariable("x", 'f', dimensions=("x",))
x_var.units = "m";
x_var[:] = x

topg_var = nc.createVariable("topg", 'f', dimensions=("x",))
topg_var.units = "m";
topg_var.standard_name = "bedrock_altitude"
topg_var[:] = topg

thk_var = nc.createVariable("thk", 'f', dimensions=("x",))
thk_var.units = "m";
thk_var.standard_name = "land_ice_thickness"
thk_var[:] = thk

qgeo = 0.042
bheatflx_var = nc.createVariable("bheatflx", 'f', dimensions=("x",))
bheatflx_var.units = "W m-2";
bheatflx_var[:] = qgeo*np.ones_like(x);

nc.close()
