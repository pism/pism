#!/usr/bin/env python

# Copyright (C) 2012, 2014 Moritz Huetten

# geometry setup MISMIP3D P75S/P10S-experiment

from pylab import *
import os
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

subgl = False

#### command line arguments ####
try:
    opts, args = getopt.getopt(sys.argv[1:], "si:a:", ["subgl", "inpath=", "amplitude="])
    for opt, arg in opts:
        if opt in ("-s", "--subgl"):
            subgl = True
        if opt in ("-i", "--inpath"):
            inpath = arg
        if opt in ("-a", "--amplitude"):
            a = arg  # pertubation amplitude for tauc


except getopt.GetoptError:
    print('Incorrect command line arguments')
    sys.exit(2)


a = float(a) * 100.0
a = int(a)

if a == 75:
    WRIT_FILE = 'MISMIP3D_P75S_initialSetup.nc'
    a = 0.75
elif a == 10:
    WRIT_FILE = 'MISMIP3D_P10S_initialSetup.nc'
    a = 0.1
else:
    WRIT_FILE = 'dummy'


######## geometry setup (moritz.huetten@pik) #########


### CONSTANTS ###

secpera = 31556926.
ice_density = 910.0

yExtent = 2 * 50  # in km
xExtent = 2 * 800  # in km


# load data from Stnd-experiment

try:
    name = inpath  # + '.nc'
    infile = NC(name, 'r')
except Exception:
    print("file '%s' not found" % name)
    sys.exit(2)
    # exit(-1)


x = squeeze(infile.variables["x"][:])
nx = len(x)


boxWidth = x[nx - 1] / ((nx - 1) / 2) / 1e3


# grid size: # of boxes

ny = int(np.floor(yExtent / boxWidth / 2) * 2 + 1)  # make it an odd number


# grid size: extent in km's, origin (0,0) in the center of the domain


y = np.linspace(-yExtent / 2, yExtent / 2, ny) * 1000.0

nxcenter = int(np.floor(0.5 * nx))
nycenter = int(np.floor(0.5 * ny))

thk = np.zeros((ny, nx))
topg = np.zeros((ny, nx))
ice_surface_temp = np.zeros((ny, nx))
precip = np.zeros((ny, nx))
tauc = np.zeros((ny, nx))

print("Informations from createSetup_PXXS.py:")
print("grid size:")
print(nx)
print(ny)
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
print("y-boxWidth:")
print(y[ny - 1] - y[ny - 2])


# load data from Stnd-result:

Thk_stnd = squeeze(infile.variables["thk"][:])
precip_stnd = squeeze(infile.variables["climatic_mass_balance"][:])
Topg_stnd = squeeze(infile.variables["topg"][:])
if subgl == True:
    Gl_mask = squeeze(infile.variables["gl_mask"][:])

print("number snapshots:")
print(len(Thk_stnd[:, 0, 0]))
lastslice = len(Thk_stnd[:, 0, 0]) - 1

thk_stnd = Thk_stnd[lastslice, :, :]
precip_stnd = precip_stnd[lastslice, :, :]
topg_stnd = Topg_stnd[lastslice, :, :]
if subgl == True:
    gl_mask = Gl_mask[lastslice, :, :]


# load bedrock geometry topg from Stnd-experiment:

for i in range(0, nx):
    for j in range(0, ny):
        topg[j, i] = topg_stnd[0, i]


# load initial ice-thickness:

for i in range(0, nx):
    for j in range(0, ny):
        thk[j, i] = thk_stnd[0, i]


# load precipitation field:

for i in range(0, nx):
    for j in range(0, ny):
        precip[j, i] = precip_stnd[0, 0] / secpera / ice_density
print("snow per year in meters")
print(precip_stnd[0, 0])


# defining dummy temperature:

for i in range(0, nx):
    for j in range(0, ny):
        ice_surface_temp[j, i] = 268.15


# number of grid cells
Mx = x.shape[0]
middle = (Mx - 1) / 2
x1 = x[middle:Mx] / 1000.0  # km
dx = x1[1] - x1[0]
thk_stnd1 = thk_stnd[1, middle:Mx]  # 1D
Mask = squeeze(infile.variables["mask"][:])
mask = Mask[lastslice, :, :]
mask = mask[1, middle:Mx]  # 1D

# find grounding line
for i in range(mask.shape[0]):
    if (thk_stnd1[i] > 0 and mask[i] == 2 and mask[i + 1] == 3):
        xg = x1[i]
        if subgl == True:
            xg_new = xg + dx / 2.0 - (1 - gl_mask[0, i]) * dx + gl_mask[0, i + 1] * dx
        else:
            xg_new = xg + dx / 2.0

print("old grounding line at position:")
print(xg, "km")

print("new grounding line at position:")
print(xg_new, "km")


xg = xg_new * 1.0e3


# defining tauc:

xb = xg
yb = 0
xc = 150e3
yc = 10e3
C = 1.0e7

a = float(a)

for i in range(nxcenter, nx):
    for j in range(0, ny):
        tauc[j, i] = C * (1 - a * math.exp(-(x[i] - xb) ** 2 / (2 * xc ** 2) - (y[j] - yb) ** 2 / (2 * yc ** 2)))

for i in range(0, nxcenter):
    for j in range(0, ny):
        tauc[j, i] = C * (1 - a * math.exp(-(x[i] + xb) ** 2 / (2 * xc ** 2) - (y[j] - yb) ** 2 / (2 * yc ** 2)))


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
        'tauc': ['Pa',
                 'yield stress for basal till (plastic or pseudo-plastic model)',
                 'yield_stress_for_basal_till',
                 1e6,
                 tauc],
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
