#!/usr/bin/env python

# Copyright (C) 2012 Ricarda Winkelmann and Torsten Albrecht and Ed Bueler

import sys
import getopt
import math
from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

WRIT_FILE = 'circular_withshelf_12km.nc'
SECPERA = 3.1556926e7


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
sizeofdomain=3600.0
dxNEW = sizeofdomain/(MxNEW-1)*1000.0 # meters
dyNEW = sizeofdomain/(MyNEW-1)*1000.0 # meters
print "dx = %.2f km, dy = %.2f km" % (dxNEW/1000.0,dyNEW/1000.0)

rho_ice	  =  910.0 # in kg/m^3
rho_ocean = 1028.0 # in kg/m^3 

accuRate= 0.3

imiddle = (MxNEW-1)/2
jmiddle = (MyNEW-1)/2

thk_max = 4000.0 # maximal ice thickness in m
thk_cf  =  200.0 # ice thickness at calving front
topg_min = -3000.0  # for practical reasons (e.g. viewing), we don't need it too deep

i_cf_m = dxNEW * (MxNEW-1) / 3 # at 1200 km
j_cf_m = dyNEW * (MyNEW-1) / 3 # at 1200 km
r_cf = math.sqrt((i_cf_m**2)+(j_cf_m**2)) # calving front position in m

gl_tol = 50.0 # tolerance for finding grounding line, depends on coarseness of the grid

# for slope of bed; shape from MISMIP
n  = 		    (729.0)
m2 = (- 2184.80/(750.0)**2)
m4 =   (1031.72/(750.0)**4)
m6 = (-  151.72/(750.0)**6)


### CREATE ARRAYS which will go in output ###
thk    = zeros((MxNEW, MyNEW), float32)		# sheet/shelf thickness
bed    = zeros((MxNEW, MyNEW), float32)		# topography -> bedrock surface elevation
Ts     = zeros((MxNEW, MyNEW), float32)		# surface temperature
accum  = zeros((MxNEW, MyNEW), float32)		# accumulation/ ablation

x=linspace(-sizeofdomain/2*1000.0,sizeofdomain/2*1000.0,MxNEW)
y=linspace(-sizeofdomain/2*1000.0,sizeofdomain/2*1000.0,MyNEW)

### BEDROCK AND ICE THICKNESS ###
for i in range(MxNEW):
  for j in range(MyNEW):
	
	### polar coordinates ###
	inew = i - imiddle
	jnew = j - jmiddle
	inew_m = dxNEW * inew # inew in meters
	jnew_m = dyNEW * jnew # jnew in meters
	radius = math.sqrt((inew_m**2)+(jnew_m**2)) # radius in m

	### set bedrock as in MISMIP experiment ###
	s = radius/1000.0
	bed[i,j] = (n + m2*s**2 + m4*s**4 + m6*s**6) # in m
	if bed[i,j] < topg_min:
	        bed[i,j] = topg_min
	
	### set thickness ###
	a = -(thk_cf - thk_max)/(r_cf)**4
	b = 2*(thk_cf - thk_max)/(r_cf)**2
	c = thk_max
	
	if (radius <= r_cf): 
		thk[i,j] = a * (radius)**4 + b* (radius)**2 + c
	else: 
		thk[i,j] = 0.0
	
	### set accumulation in m/s ###
	accum[i,j] = accuRate / SECPERA
	
	### set surface temperature ###
	Ts[i,j] = 247.0


### GROUNDING LINE RADIUS ###
for i in range(MxNEW):
	# flotation criterion: grounding line is where
	#    bed(r_gl) + thk(r_gl) = (1.0 - rho_ice/rho_ocean) * thk (r_gl)
	hgrounded = bed[i,jmiddle] + thk[i,jmiddle]
	hfloating = (1.0 - rho_ice / rho_ocean) * thk [i,jmiddle]	
	if (abs(hgrounded - hfloating) < gl_tol):
		#print "  grounding line found at ", i, jmiddle
		r_gl = math.sqrt((dxNEW*(i - imiddle))**2) # radius in m
print "grounding line radius = %.2f km" % (r_gl/1000.0)


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
	'ice_surface_temp': ['K',
             'annual mean air temperature at ice surface',
	     'surface_temperature',
	      248.0,
	      Ts],
	'climatic_mass_balance': ['m s-1',
	     'mean annual net ice equivalent accumulation rate',
             'land_ice_surface_specific_mass_balance',
	      0.2/SECPERA,
	      accum],
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

