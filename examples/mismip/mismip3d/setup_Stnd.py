#!/usr/bin/env python

# Copyright (C) 2012, 2014 Moritz Huetten

import sys
import getopt
import math
#from numpy import *
import numpy as np
try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


# geometry setup MISMIP3D Stnd-experiment


WRIT_FILE = 'MISMIP3D_Stnd_initialSetup.nc'

accumrate = 0.5  # accumulation rate in m/a

#### command line arguments ####
try:
    opts, args = getopt.getopt(sys.argv[1:], "a:r:", ["accumrate=", "resolution="])
    for opt, arg in opts:
        if opt in ("-a", "--accumulation"):
            accumrate = arg
        if opt in ("-r", "--resolution"):  # resolution in km
            boxWidth = arg
except getopt.GetoptError:
    print('Incorrect command line arguments')
    sys.exit(2)


######## geometry setup (moritz.huetten@pik) #########


boxWidth = float(boxWidth)
accumrate = float(accumrate)

### CONSTANTS ###


secpera = 31556926.
ice_density = 910.0             # [kg m-3]

yExtent = 2 * boxWidth  # in km
xExtent = 2 * 800  # in km


# grid size: # of boxes

ny = int(np.floor(yExtent / boxWidth / 2) * 2 + 1)  # make it an odd number
nx = int(np.floor(xExtent / boxWidth / 2) * 2 + 1)  # make it an odd number

# grid size: extent in km's, origin (0,0) in the center of the domain

x = np.linspace(-xExtent / 2, xExtent / 2, nx) * 1000.0
y = np.linspace(-yExtent / 2, yExtent / 2, ny) * 1000.0

nxcenter = int(np.floor(0.5 * nx))
nycenter = int(np.floor(0.5 * ny))

thk = np.zeros((ny, nx))
topg = np.zeros((ny, nx))
ice_surface_temp = np.zeros((ny, nx))
precip = np.zeros((ny, nx))

# define bed elevation:
# linear sloping bed in x-direction with top in the middle, symmetric in y-direction

# print range(0,int(np.floor(0.5*nx)))

print("Informations from createSetup_Stnd.py:")
print("grid size (pixel):")
print(ny)
print(nx)
print("grid size center:")
print(nxcenter)
print(nycenter)

print("domain range in meters:")
print("x-dir:")
print(x[0])
print(x[nx - 1])
print("y-dir:")
print(y[0])
print(y[ny - 1])


# define bedrock geometry topg:

for i in range(0, nx):
    for j in range(0, ny):
        topg[j, i] = -100.0 - abs(x[i]) / 1.0e3

# define constant initial ice-thickness and extent:

thickness = 500.0  # initial ice thickness in meters
xfront = 700.0  # x-position of fixed calving front in km

nxfront = int(xfront / boxWidth)

for i in range(nxcenter - nxfront, nxcenter + nxfront):
    for j in range(0, ny):
        thk[j, i] = thickness


# define precipitation field:

for i in range(0, nx):
    for j in range(0, ny):
        precip[j, i] = accumrate / secpera


# defining dummy temperature:

for i in range(0, nx):
    for j in range(0, ny):
        ice_surface_temp[j, i] = 268.15


##### define dimensions in NetCDF file #####
ncfile = NC(WRIT_FILE, 'w', format='NETCDF3_CLASSIC')
xdim = ncfile.createDimension('x', nx)
ydim = ncfile.createDimension('y', ny)

##### define variables, set attributes, write data #####
# format: ['units', 'long_name', 'standard_name', '_FillValue', array]

vars = {'y':   	['m',
                 'y-coordinate in Cartesian system',
                 'projection_y_coordinate',
                 None,
                 y],
        'x':   	['m',
                 'x-coordinate in Cartesian system',
                 'projection_x_coordinate',
                 None,
                 x],
        'thk': 	['m',
                 'floating ice shelf thickness',
                 'land_ice_thickness',
                 1.0,
                 thk],
        'topg': ['m',
                 'bedrock surface elevation',
                 'bedrock_altitude',
                 -600.0,
                 topg],
        'ice_surface_temp': ['K',
                             'annual mean air temperature at ice surface',
                             'surface_temperature',
                             248.0,
                             ice_surface_temp],
        'climatic_mass_balance': ['kg m-2 year-1',
                                  'mean annual net ice equivalent accumulation rate',
                                  'land_ice_surface_specific_mass_balance_flux',
                                  0.2 * ice_density,
                                  precip],
        }

for name in list(vars.keys()):
    [_, _, _, fill_value, data] = vars[name]
    if name in ['x', 'y']:
        var = ncfile.createVariable(name, 'f4', (name,))
    else:
        var = ncfile.createVariable(name, 'f4', ('y', 'x'), fill_value=fill_value)
    for each in zip(['units', 'long_name', 'standard_name'], vars[name]):
        if each[1]:
            setattr(var, each[0], each[1])
    var[:] = data

# finish up
ncfile.close()
print("NetCDF file ", WRIT_FILE, " created")
print('')
