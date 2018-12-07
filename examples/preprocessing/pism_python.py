#!/usr/bin/env python

# Copyright (C) 2009-2015, 2018 the PISM Authors

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

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
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
acab = np.zeros((Mx, My))
artm = np.zeros((Mx, My)) + 273.15 + 10.0  # 10 degrees Celsius
topg = 1000.0 + 200.0 * (xx + yy) / max(Lx, Ly)  # change "1000.0" to "0.0" to test
# flotation criterion, etc.
thk = 3000.0 * (1.0 - 3.0 * (xx ** 2 + yy ** 2) / Lx ** 2)
thk[thk < 0.0] = 0.0

# Output filename
ncfile = 'foo.nc'

# Write the data:
nc = CDF(ncfile, "w", format='NETCDF3_CLASSIC')  # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=Mx)
nc.createDimension("y", size=My)

x_var = nc.createVariable("x", 'f4', dimensions=("x",))
x_var.units = "m"
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = x

y_var = nc.createVariable("y", 'f4', dimensions=("y",))
y_var.units = "m"
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"
y_var[:] = y

fill_value = np.nan


def def_var(nc, name, units, fillvalue):
    # dimension transpose is standard: "float thk(y, x)" in NetCDF file
    var = nc.createVariable(name, 'f', dimensions=("y", "x"), fill_value=fillvalue)
    var.units = units
    return var


bed_var = def_var(nc, "topg", "m", fill_value)
bed_var.standard_name = "bedrock_altitude"
bed_var[:] = topg

thk_var = def_var(nc, "thk", "m", fill_value)
thk_var.standard_name = "land_ice_thickness"
thk_var[:] = thk

acab_var = def_var(nc, "climatic_mass_balance", "m year-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = acab

artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
artm_var[:] = artm

# set global attributes
nc.Conventions = "CF-1.4"
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()
print('  PISM-bootable NetCDF file %s written' % ncfile)
print('  for example, run:')
print('    $ pismr -i foo.nc -bootstrap -Mx 41 -My 41 -Mz 21 -Lz 4000 -Mbz 5 -Lbz 500 -y 1')
