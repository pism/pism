#!/usr/bin/env python

# Copyright (C) 2013 the PISM Authors

# This script sets up the bootstrap file.  It also converts the configuration
# parameter file into NetCDF:
#   gumparams.cdl  -->  gumparams.nc
# See also the run script:
#   $ ./rungum.sh

import sys
import time
import numpy as np

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

# lab setup is table with hole in the middle into which is piped the
# shear-thinning fluid, which is Xanthan gum 1% solution
Lx = 250.0e-3    # m;  = 250 mm;  table is approx 500 mm x 500 mm?
Ly = Lx          # square table
flux = 3.0e-3    # kg s-1;  = 3 g s-1;  
pipeR = 20.0e-3  # m;  = 20 mm;  input pipe has this radius  FIXME: GUESS
rho = 1000.0     # kg m-3;  density of gum = density of fresh water
temp = 20.0      # C;  fluid is at 20 deg

# see gumparams.cdl for additional parameter settings
from subprocess import call
CONF = 'gumparams'
call(['rm', '-f', CONF + '.nc'])
call(['ncgen', '-o', CONF + '.nc', CONF + '.cdl'])
print('  PISM-readable config override file %s written' % (CONF + '.nc'))

# set up the grid:
Mx = 500
My = 500
x = np.linspace(-Lx,Lx,Mx)
y = np.linspace(-Ly,Ly,My)

# create dummy fields
[xx,yy] = np.meshgrid(x,y);  # if there were "ndgrid" in numpy we would use it
topg = np.zeros((Mx,My))
thk  = np.zeros((Mx,My))  # no fluid on table at start
artm = np.zeros((Mx,My)) + 273.15 + temp; # 20 degrees Celsius
fluxthickness = flux / (rho * np.pi * pipeR**2)  # flux as m s-1
acab = np.zeros((Mx,My)) + fluxthickness;
acab[xx**2 + yy**2 > pipeR**2] = 0.0;

# Output filename
ncfile = 'initgum.nc'

# Write the data:
nc = CDF(ncfile, "w",format='NETCDF3_CLASSIC') # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=Mx)
nc.createDimension("y", size=My)

x_var = nc.createVariable("x", 'f4', dimensions=("x",))
x_var.units = "m";
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = x

y_var = nc.createVariable("y", 'f4', dimensions=("y",))
y_var.units = "m";
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

acab_var = def_var(nc, "climatic_mass_balance", "m s-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = acab

artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
artm_var[:] = artm

# set global attributes
nc.Conventions = 'CF-1.4'
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()
print('  PISM-bootable NetCDF file %s written' % ncfile)
print('  now run "rungum.sh"')
print('  TRY:  $ pismr -config_override gumparams.nc -boot_file initgum.nc -no_energy -cold -sia_flow_law isothermal_glen -sia_e 1.0 -Mx 101 -My 101 -Mz 51 -Lz 0.050 -Mbz 0 -Lbz 0 -z_spacing equal -y 1e-7 -max_dt 1e-9')

