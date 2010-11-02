#!/usr/bin/env python

# Copyright (C) 2009-2010 Andy Aschwanden
# Author: 	Andy Aschwanden
#		Arctic Region Supercomputing Center
#		University of Alaska Fairbanks
#		andy.aschwanden@arsc.edu

## @package pism_python
# \author Andy Aschwanden, University of Alaska Fairbanks, USA
# \brief Creates "from scratch" a very boring dataset with the right format
# to use as a PISM bootstrapping file.
# \details Illustration use of Python for this purpose.
#
# Usage, including the minimal kind of PISM call needed to bootstrap from
# this file:
#
# \verbatim $ pism_python.py  # creates foo.nc \endverbatim
# \verbatim $ pismr -boot_file foo.nc -Mx 101 -My 201 -surface constant \
#                   -Mz 11 -Lz 4000 -Mbz 3 -Lbz 2000 -y 1  \endverbatim

import sys
import time
import numpy as np

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

# set up the grid:
Lx = 1e6;
Ly = 1e6;
Mx = 101;
My = 201;
x = np.linspace(-Lx,Lx,Mx);
y = np.linspace(-Ly,Ly,My);

# create dummy fields
[yy,xx] = np.meshgrid(x,y);
topg = (xx + yy) / max(Lx, Ly);
acab = np.zeros((My,Mx));
artm = np.zeros((My,Mx)) + 10.0; # 10 degrees Celsius
thk  = np.zeros((My,Mx)) + 1.0; # 1 km thick

# Output filename
ncfile = 'foo.nc'

# Write the data:
nc = CDF(ncfile, "w",format='NETCDF3_CLASSIC') # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=Mx)
nc.createDimension("y", size=My)

x_var = nc.createVariable("x", 'f4', dimensions=("x",))
x_var.units = "m";
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"

y_var = nc.createVariable("y", 'f4', dimensions=("y",))
y_var.units = "m";
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"

x_var[:] = x
y_var[:] = y

fill_value = np.nan

def def_var(nc, name, units, fillvalue):
    var = nc.createVariable(name, 'f', dimensions=("y", "x"), fill_value=fillvalue)
    var.units = units
    return var


bed_var = def_var(nc, "topg", "m", fill_value)
bed_var.standard_name = "bedrock_altitude"
bed_var[:] = topg

thk_var = def_var(nc, "thk", "m", fill_value)
thk_var.standard_name = "land_ice_thickness"
thk_var[:] = np.flipud(thk)

acab_var = def_var(nc, "acab", "m year-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = np.flipud(acab)

artm_var = def_var(nc, "artm", "K", fill_value)
artm_var[:] = np.flipud(artm)

# set global attributes
nc.Conventions = "CF-1.4"
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()
print('Done writing NetCDF file %s!\n' % ncfile)        
