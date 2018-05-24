#!/usr/bin/env python
# Copyright (C) 2011, 2013, 2014, 2016, 2018 Torsten Albrecht and Moritz Huetten

# ./createSetup_flowline.py -a 0.0 -r 10.0

import sys
import getopt
import math
import numpy as np
try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

# geometry setup flowline

WRIT_FILE = 'flowline_setup.nc'  # default name
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


# geometry setup
boxWidth = float(boxWidth)
accumrate = float(accumrate)
WRIT_FILE = 'flowline_setup_' + str(int(boxWidth)) + 'km.nc'

### CONSTANTS ###
secpera = 31556926.
yExtent = 2 * boxWidth  # in km
xExtent = 500  # in km
standard_gravity = 9.81                 # g
B0 = 1.9e8                              # ice hardness
rho_ice = 910.0  # in kg/m^3
rho_ocean = 1028.0  # in kg/m^3

# create config overwrite file as used in T. Albrecht, M. A. Martin, R. Winkelmann, M. Haseloff, A. Levermann; Parameterization for subgrid-scale motion of ice-shelf calving fronts; (2011), The Cryosphere 5, 35-44, DOI:10.5194/tc-5-35-2011.
'''Generates a config file containing flags and parameters
for a particular experiment and step.

This takes care of flags and parameters that *cannot* be set using
command-line options. (We try to use command-line options whenever we can.)
'''

filename = "flowline_config.nc"
nc = NC(filename, 'w', format="NETCDF3_CLASSIC")
var = nc.createVariable("pism_overrides", 'i')
attrs = {"ocean.always_grounded": "no",
         "geometry.update.use_basal_melt_rate": "no",
         "stress_balance.ssa.compute_surface_gradient_inward": "no",
         "flow_law.isothermal_Glen.ice_softness": (B0) ** -3,
         "constants.ice.density": rho_ice,
         "constants.sea_water.density": rho_ocean,
         "bootstrapping.defaults.geothermal_flux": 0.0,
         "stress_balance.ssa.Glen_exponent": 3.0,
         "constants.standard_gravity": standard_gravity,
         "ocean.sub_shelf_heat_flux_into_ice": 0.0,
         "stress_balance.sia.bed_smoother.range": 0.0,
         }

for name, value in attrs.items():
    var.setncattr(name, value)
nc.close()


# grid size: # of boxes

ny = int(np.floor(yExtent / boxWidth / 2) * 2 + 1)  # make it an odd number
nx = int(np.floor(xExtent / boxWidth / 2) * 2 + 1)  # make it an odd number

# grid size: extent in km's, origin (0,0) to the left of the domain

x = np.linspace(0, xExtent, nx) * 1000.0
y = np.linspace(-yExtent / 2, yExtent / 2, ny) * 1000.0

nxcenter = int(np.floor(0.5 * nx))
nycenter = int(np.floor(0.5 * ny))

thk = np.zeros((ny, nx))
topg = np.zeros((ny, nx))
ice_surface_temp = np.zeros((ny, nx))
precip = np.zeros((ny, nx))


print("Informations from createSetup_flowline.py:")
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


thickness = 600.0  # initial ice thickness in meters
xfront = 700.0  # x-position of fixed calving front in km

# bc
distbc = 50.0  # km
igl = int(np.floor(distbc / boxWidth))
vel_bc = 300  # m/yr
bc_mask = np.zeros((ny, nx))
ubar = np.zeros((ny, nx))
vbar = np.zeros((ny, nx))

# "typical constant ice parameter" as defined in the paper and in Van der
# Veen's "Fundamentals of Glacier Dynamics", 1999
#C = (rho_ice * standard_gravity * (1.0 - rho_ice/rho_ocean) / (4 * B0))**3
# H0=thickness
# Q0=vel_bc*thickness/secpera

nxfront = int(xfront / boxWidth)

# define bedrock geometry topg:

for i in range(0, nx):
    for j in range(0, ny):
        topg[j, i] = -2000.0

# define constant initial ice-thickness and extent:

# for i in range(nxcenter-nxfront,nxcenter+nxfront):
for i in range(0, nx):
    for j in range(0, ny):
        if i <= igl:
            #thk[j,i] = (4.0 * C / Q0 * (i-igl) + 1 / H0**4)**(-0.25)
            thk[j, i] = thickness

# define precipitation field:

for i in range(0, nx):
    for j in range(0, ny):
        precip[j, i] = accumrate / secpera

# defining dummy temperature:

for i in range(0, nx):
    for j in range(0, ny):
        ice_surface_temp[j, i] = 247.0

# defining boundary conditions:

for i in range(0, nx):
    for j in range(0, ny):
        if i <= igl:
            ubar[j, i] = vel_bc / secpera
            vbar[j, i] = 0.0
            bc_mask[j, i] = 1.0


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
                                  0.2 * 910.0,
                                  precip],
        'bc_mask': ['',
                    'bc_mask',
                    'bc_mask',
                    0.0,
                    bc_mask],
        'u_ssa_bc': ['m s-1',
                     'X-component of the SSA velocity boundary conditions',
                     'ubar',
                     0.0,
                     ubar],
        'v_ssa_bc': ['m s-1',
                     'Y-component of the SSA velocity boundary conditions',
                     'vbar',
                     0.0,
                     vbar],
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
