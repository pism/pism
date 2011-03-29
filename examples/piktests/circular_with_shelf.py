#!/usr/bin/env python

# Copyright (C) 2011 Ricarda Winkelmann and Torsten Albrecht

import sys
import getopt
import math
import time
from numpy import *
from netCDF4 import Dataset as NC


#######################
#### set constants ####
#######################
WRIT_FILE = 'circular_with_shelf_12km.nc'
VERBOSE = 0
SECPERA = 3.1556926e7

################################################################
#### function which will print ignored lines if VERBOSE > 0 ####
################################################################
def vprint(s):
  if VERBOSE > 0:
    print s

######################################################################################
#### function to convert a (Mx,My) array into a (1,Mx,My) array; for NetCDF write ####
######################################################################################
def addTime(A):
  shA = shape(A)
  shB = (1,shA[0],shA[1])
  return reshape(A,shB)

################################
#### command line arguments ####
################################
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:o:v:", ["prefix=", "out=", "verbose="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      BED_FILE = arg + DATA_FILE
    if opt in ("-o", "--out"):
      WRIT_FILE = arg
    if opt in ("-v", "--verbose"):
      verbose = float(arg)
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)



################################################
######## geometry setup (ricardaw@pik) #########
################################################

### CONSTANTS ###
MxNEW = 301 # should be odd number
MyNEW = 301
sizeofdomain=3600.0
dxNEW = sizeofdomain/(MxNEW-1)*1000.0 # meters
dyNEW = sizeofdomain/(MyNEW-1)*1000.0 # meters
oceanborder = 3

rho_ice	  =  910.0 # in kg/m^3
rho_ocean = 1028.0 # in kg/m^3 

accuRate= 0.3
myshift = 0.0
myscale = 1.0

imiddle = (MxNEW-1)/2
jmiddle = (MyNEW-1)/2

thk_max = 4000.0 # maximal ice thickness in m
thk_cf  =  200.0 # ice thickness at calving front

i_cf_m = dxNEW * (MxNEW-1) / 6 # at 600 km
j_cf_m = dyNEW * (MyNEW-1) / 6 # at 600 km
r_cf = math.sqrt((i_cf_m**2)+(j_cf_m**2)) # calving front position in m

gl_tol = 50.0 # tolerance for finding grounding line, depends on coarseness of the grid

# for slope of bed
n  = 		    (729.0)
m2 = (- 2184.80/(750.0)**2)
m4 =   (1031.72/(750.0)**4)
m6 = (-  151.72/(750.0)**6)


### CREATE ARRAYS ###
radius = zeros((MxNEW, MyNEW), float32)		# sheet / shelf is going to have circular symmetry
hgrounded = zeros((MxNEW, MyNEW), float32)
hfloating = zeros((MxNEW, MyNEW), float32)
mask   = zeros((MxNEW,MyNEW), int16)		# 1->grounded, 2->dragging, 3->floating, 7->floating_ocean0
thk    = zeros((MxNEW, MyNEW), float32)		# sheet/shelf thickness
bed    = zeros((MxNEW, MyNEW), float32)		# topography -> bedrock surface elevation
Ts     = zeros((MxNEW, MyNEW), float32)		# surface temperature
accum  = zeros((MxNEW, MyNEW), float32)		# accumulation/ ablation

x=linspace(-sizeofdomain/2*1000.0,sizeofdomain/2*1000.0,MxNEW)
y=linspace(-sizeofdomain/2*1000.0,sizeofdomain/2*1000.0,MyNEW)

### BEDROCK AND ICE THICKNESS ###
for i in range(MxNEW):
  for j in range(MyNEW):
	
	### set border of computational domain to ocean always ###
	if (i >= MxNEW-oceanborder): mask[i,j] = 7
	if (i < oceanborder): mask[i,j] = 7
	if (j >= MyNEW-oceanborder): mask[i,j] = 7
	if (j < oceanborder): mask[i,j] = 7
	
	### polar coordinates ###
	inew = i - imiddle
	jnew = j - jmiddle
	inew_m = dxNEW * inew # inew in meters
	jnew_m = dyNEW * jnew # jnew in meters
	radius[i,j] = math.sqrt((inew_m**2)+(jnew_m**2)) # radius in m

	### set bedrock as in MISMIP experiment ###
	bed[i,j] = (n + m2*(radius[i,j]/1000.0)**2 + m4*(radius[i,j]/1000.0)**4 + m6*(radius[i,j]/1000.0)**6) # in m
	
	### set thickness ###
	a = -(thk_cf - thk_max)/(r_cf)**4
	b = 2*(thk_cf - thk_max)/(r_cf)**2
	c = thk_max
	
	if (radius[i,j] <= r_cf): 
		thk[i,j] = a * (radius[i,j])**4 + b* (radius[i,j])**2 + c
	else: 
		thk[i,j] = 0.0
	
	### set accumulation in m/s ###
	accum[i,j] = accuRate / (365 * 24 * 60 * 60)
	
	### set surface temperature ###
	Ts[i,j] = (247 + myshift) * myscale


### GROUNDING LINE RADIUS ###
for i in range(MxNEW):
	# floatation criterion: grounding line is where bed(r_gl) + thk(r_gl) = (1.0 - rho_ice/rho_ocean) * thk (r_gl)
	hgrounded[i,jmiddle] = bed[i,jmiddle] + thk[i,jmiddle]
	hfloating[i,jmiddle] = (1.0 - rho_ice / rho_ocean) * thk [i,jmiddle]	
	if (abs(hgrounded[i,jmiddle] - hfloating[i,jmiddle]) < gl_tol):
		print "habs bei ", i, jmiddle
		r_gl = math.sqrt((dxNEW*(i - imiddle))**2) # radius in m	
	

### MASK ###
for i in range(MxNEW):
  for j in range(MyNEW):
	if (radius[i,j] <= r_gl): mask[i,j] = 2
	if ((radius[i,j] > r_gl) and radius[i,j] <= r_cf): mask[i,j] = 3
	



#####################################################################
##### create and define dimensions and variables in NetCDF file #####
#####################################################################

##### define dimensions in NetCDF file #####
ncfile = NC(WRIT_FILE, 'w',format='NETCDF3_CLASSIC')
xdim = ncfile.createDimension('x', MyNEW)
ydim = ncfile.createDimension('y', MxNEW)



##### define variables, set attributes, write data #####
# format: ['units', 'long_name', 'standard_name', '_FillValue', array]
vars = {'y': ['m',
              'x-coordinate in Cartesian system',
              'projection_x_coordinate',
              None,
              x],
        'x': ['m',
              'y-coordinate in Cartesian system',
              'projection_y_coordinate',
              None,
              y],
		'thk': ['m',
		        'floating ice shelf thickness',
		        'land_ice_thickness',
		         1.0,
		         thk],
		'radius': ['m',
				 'radius',
				 'radius',
				  1.0,
				 radius],
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
        'mask': [None,
                 'grounded or floating integer mask',
                 None,
                 None,
                 mask],}


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

##### attributes in NetCDF file #####
# set global attributes
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(ncfile, 'history', historystr)
setattr(ncfile, 'Conventions', 'CF-1.4') # only global attribute

# finish up
ncfile.close()
print "NetCDF file ",WRIT_FILE," created"

