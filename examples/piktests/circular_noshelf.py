#!/usr/bin/env python

# Copyright (C) 2011 Ricarda Winkelmann and Torsten Albrecht and Ed Bueler

import sys
import getopt
import math
from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

WRIT_FILE = 'circular_noshelf.nc'
SECPERA = 3.15569259747e7               # seconds per yea
standard_gravity = 9.81                 # g
B0 = 1.9e8                              # ice hardness

rho_ice	  =  910.0 # in kg/m^3
rho_ocean = 1028.0 # in kg/m^3
# "typical constant ice parameter" as defined in the paper and in Van der
# Veen's "Fundamentals of Glacier Dynamics", 1999
C = (rho_ice * standard_gravity * (1.0 - rho_ice/rho_ocean) / (4 * B0))**3


#### command line arguments ####
try:
  opts, args = getopt.getopt(sys.argv[1:], "o:", ["out="])
  for opt, arg in opts:
    if opt in ("-o", "--out"):
      WRIT_FILE = arg
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)


######## geometry setup (ricardaw@pik) #########

### CONSTANTS ###
MxNEW = 301 # should be odd number
MyNEW = 301
sizeofdomain=1000.0
dxNEW = sizeofdomain/(MxNEW-1)*1000.0 # meters
dyNEW = sizeofdomain/(MyNEW-1)*1000.0 # meters
print "dx = %.2f km, dy = %.2f km" % (dxNEW/1000.0,dyNEW/1000.0)

vel_bc=300 # m/yr
accuRate= 0.3

imiddle = (MxNEW-1)/2
jmiddle = (MyNEW-1)/2

topg_min = -3000.0  # for practical reasons (e.g. viewing), we don't need it too deep

r_cf = 0.3e6
gl_tol = 50.0 # tolerance for finding grounding line, depends on coarseness of the grid



### CREATE ARRAYS which will go in output ###
thk    = zeros((MxNEW, MyNEW), float32)		# sheet/shelf thickness
bed    = zeros((MxNEW, MyNEW), float32)		# topography -> bedrock surface elevation
Ts     = zeros((MxNEW, MyNEW), float32)		# surface temperature
accum  = zeros((MxNEW, MyNEW), float32)		# accumulation/ ablation

bcflag  = zeros((MxNEW, MyNEW), float32)		# accumulation/ ablation
ubar  = zeros((MxNEW, MyNEW), float32)		# accumulation/ ablation
vbar  = zeros((MxNEW, MyNEW), float32)		# accumulation/ ablation

x=linspace(-sizeofdomain/2*1000.0,sizeofdomain/2*1000.0,MxNEW)
y=linspace(-sizeofdomain/2*1000.0,sizeofdomain/2*1000.0,MyNEW)

### GROUNDING LINE RADIUS ###
r_gl = 0.25e6


### BEDROCK AND ICE THICKNESS ###
for i in range(MxNEW):
  for j in range(MyNEW):
	
	### polar coordinates ###
	inew = i - imiddle
	jnew = j - jmiddle
	inew_m = dxNEW * inew # inew in meters
	jnew_m = dyNEW * jnew # jnew in meters
	radius = math.sqrt((inew_m**2)+(jnew_m**2)) # radius in m

	if radius < r_gl:	
		bed[i,j]=100.0
	else: 
		bed[i,j]=topg_min
	
	H0=600
	Q0=300*600/SECPERA

	
	if (radius <= r_cf and radius >r_gl): 
		#thk[i,j] = (4.0 * C / Q0 * (radius-r_gl) + 1 / H0**4)**(-0.25)
		if (thk[i,j]>600.0):
			thk[i,j]=600.0
	elif (radius <= r_gl):
		thk[i,j] = 0.0
	else: 
		thk[i,j] = 0.0
	
	### set accumulation in m/s ###
	if (radius > r_gl):
		accum[i,j] = accuRate / SECPERA
	
	### set surface temperature ###
	Ts[i,j] = 247.0


### set bcmask and dirichelt boundary conditions ###
for i in range(MxNEW):
  for j in range(MyNEW):
	### polar coordinates ###
	inew = i - imiddle
	jnew = j - jmiddle
	inew_m = dxNEW * inew # inew in meters
	jnew_m = dyNEW * jnew # jnew in meters
	radius = math.sqrt((inew_m**2)+(jnew_m**2)) # radius in m
	width = dxNEW*3
	
	if (radius <= r_gl-width): 
		bcflag[i,j] = 1.0
	elif (radius <= r_gl):
		bcflag[i,j] = 1.0
		ubar[i,j] = vel_bc/SECPERA*inew_m/radius
		vbar[i,j] = vel_bc/SECPERA*jnew_m/radius
		thk[i,j] = 600.0
		bed[i,j] = topg_min

	else: 
		bcflag[i,j] = 0.0


##### define dimensions in NetCDF file #####
ncfile = NC(WRIT_FILE, 'w',format='NETCDF3_CLASSIC')
xdim = ncfile.createDimension('x', MyNEW)
ydim = ncfile.createDimension('y', MxNEW)

##### define variables, set attributes, write data #####
# format: ['units', 'long_name', 'standard_name', '_FillValue', array]
vars = {'x': ['m',
              'x-coordinate in Cartesian system',
              'projection_x_coordinate',
              None,
              x],
        'y': ['m',
              'y-coordinate in Cartesian system',
              'projection_y_coordinate',
              None,
              y],
	'thk': ['m',
	      'floating ice shelf thickness',
	      'land_ice_thickness',
	      1.0,
	      thk],
	'topg': ['m',
	      'bedrock surface elevation',
	      'bedrock_altitude',
	      -600.0,
	      bed],
     	'artm': ['K',
	      'annual mean air temperature at ice surface',
	      'surface_temperature',
	      248.0,
	      Ts],
	'acab': ['m s-1',
	      'mean annual net ice equivalent accumulation rate',
	      'land_ice_surface_specific_mass_balance',
	      0.2/SECPERA,
	      accum],
	'bcflag': ['',
		  'bcflag',
		  'bcflag',
		   0.0,
		   bcflag],
    'ubar': ['m s-1',
	  	  'ubar',
		  'ubar',
		   0.0,
		   vbar],
	'vbar': ['m s-1',
		   'vbar',
		   'vbar',
		    0.0,
		    ubar],
	}


for name in vars.keys():
    [_, _, _, fill_value, data] = vars[name]
    if name in ['x', 'y']:
        var = ncfile.createVariable(name, 'f4', (name,))
    else:
        var = ncfile.createVariable(name, 'f4', ('y', 'x'), fill_value = fill_value)
    for each in zip(['units', 'long_name', 'standard_name'], vars[name]):
        if each[1]:
            setattr(var, each[0], each[1])
    var[:] = data

# finish up
ncfile.close()
print "NetCDF file ",WRIT_FILE," created"

