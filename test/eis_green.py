#! /usr/bin/env python

#import Numeric
import sys
import getopt
from numpy import *
from pycdf import *

SECPERA = 3.1556926e7
GRID_FILE = 'grid20-EISMINT'
SUAQ_FILE = 'suaq20-EISMINT'
WRIT_FILE = 'eis_green20.nc'

##### command line arguments #####

try:
  opts, args = getopt.getopt(sys.argv[1:], "g:", ["grid"])
  for opt, arg in opts:
    if opt in ("-g", "--grid"):
      GRID_FILE = "grid" + opt + "-EISMINT"
      SUAQ_FILE = "suaq" + opt + "-EISMINT"
      WRIT_FILE = "eis_green" + opt + ".nc"
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

##### grid20-EISMINT #####

# open the data file and begin reading it
input=open(GRID_FILE, 'r')

dim=[]
for num in input.readline().split():
        dim.append(num)
input.readline() # get rid of the titles

lat = zeros((int(dim[0]), int(dim[1])), float32)
lon = zeros((int(dim[0]), int(dim[1])), float32)

latcount=0;
loncount=0
x=0
# read in all the data
for line in input.read().split():
  if x%4 == 2: # surface elevation
    lat[latcount%int(dim[0]), latcount/int(dim[0])]=float(line)
    latcount = latcount + 1
  elif x%4 == 3: # thickness
    lon[loncount%int(dim[0]), loncount/int(dim[0])]=float(line)
    loncount = loncount + 1
  x = x+1

# done reading from data file
input.close()

##### suaq20-EISMINT #####

# open the data file and begin reading it
input=open(SUAQ_FILE, 'r')

dim=[]
for num in input.readline().split():
	dim.append(num)
input.readline() # get rid of the titles

S = zeros((int(dim[0]), int(dim[1])), float32)
H = zeros((int(dim[0]), int(dim[1])), float32)
B = zeros((int(dim[0]), int(dim[1])), float32)
acc = zeros((int(dim[0]), int(dim[1])), float32)

Scount=0;
Hcount=0
Bcount=0
Acccount=0
x=0
# read in all the data
for line in input.read().split():
	if x%4 == 0: # surface elevation
		S[Scount%int(dim[0]), Scount/int(dim[0])]=float(line)
		Scount = Scount + 1
	elif x%4 == 1: # thickness
		H[Hcount%int(dim[0]), Hcount/int(dim[0])]=float(line)
		Hcount = Hcount + 1
	elif x%4 == 2: # bedrock
		B[Bcount%int(dim[0]), Bcount/int(dim[0])]=float(line)
		Bcount = Bcount + 1
	else: # accumulation (m/a -> m/s)
		acc[Acccount%int(dim[0]), Acccount/int(dim[0])]=float(line)/SECPERA
		Acccount = Acccount + 1
	x = x+1

# done reading from data file
input.close()

print "Total Values Read: "+str(x)

# ready to write NetCDF file
ncfile = CDF(WRIT_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()

# define the dimensions
xdim = ncfile.def_dim('x', int(dim[0]))
ydim = ncfile.def_dim('y', int(dim[1]))
tdim = ncfile.def_dim('t', NC.UNLIMITED)

# define the variables
polarVar = ncfile.def_var('polar_stereographic', NC.INT)
xvar = ncfile.def_var('x', NC.FLOAT, (xdim,))
yvar = ncfile.def_var('y', NC.FLOAT, (ydim,))
tvar = ncfile.def_var('t', NC.FLOAT, (tdim,))
lonvar = ncfile.def_var('lon', NC.FLOAT, (xdim, ydim))
latvar = ncfile.def_var('lat', NC.FLOAT, (xdim, ydim))
hvar = ncfile.def_var('h', NC.FLOAT, (xdim, ydim))
Hvar = ncfile.def_var('H', NC.FLOAT, (xdim, ydim))
Bvar = ncfile.def_var('b', NC.FLOAT, (xdim, ydim)) 
Accvar = ncfile.def_var('accum', NC.FLOAT, (xdim, ydim))

# set global attributes
setattr(ncfile, 'Conventions', 'CF-1.0')

# set the attributes of the variables
setattr(polarVar, 'grid_mapping_name', 'polar_stereographic')
setattr(polarVar, 'straight_vertical_longitude_from_pole', -40)
setattr(polarVar, 'latitude_of_projection_origin', 75)
setattr(polarVar, 'standard_parallel', 90)

setattr(xvar, 'axis', 'X')
setattr(xvar, 'long_name', 'x-coordinate in Cartesian system')
setattr(xvar, 'standard_name', 'projection_x_coordinate')
setattr(xvar, 'units', 'm')

setattr(yvar, 'axis', 'Y')
setattr(yvar, 'long_name', 'y-coordinate in Cartesian system')
setattr(yvar, 'standard_name', 'projection_y_coordinate')
setattr(yvar, 'units', 'm')

setattr(tvar, 'units', 'seconds since 2007-01-01 00:00:00')

setattr(lonvar, 'long_name', 'longitude')
setattr(lonvar, 'standard_name', 'longitude')
setattr(lonvar, 'units', 'degrees_east')

setattr(latvar, 'long_name', 'latitude')
setattr(latvar, 'standard_name', 'latitude')
setattr(latvar, 'units', 'degrees_north')

setattr(hvar, 'long_name', 'surface_altitude')
setattr(hvar, 'standard_name', 'surface_altitude')
setattr(hvar, 'units', 'm')

setattr(Hvar, 'long_name', 'land_ice_thickness')
setattr(Hvar, 'standard_name', 'land_ice_thickness')
setattr(Hvar, 'units', 'm')

setattr(Bvar, 'long_name', 'bedrock_altitude')
setattr(Bvar, 'standard_name', 'bedrock_altitude')
setattr(Bvar, 'units', 'm')
setattr(Bvar, 'missing_value', 0.0)

setattr(Accvar, 'long_name', 'mean ice equivalent accumulation rate')
setattr(Accvar, 'standard_name', 'land_ice_surface_specific_mass_balance')
setattr(Accvar, 'units', 'm s-1')

# write the data to the NetCDF file
spacing = float(dim[2])*1000
for i in range(int(dim[0])):
	xvar[i]=((1-float(dim[0]))/2+i)*spacing
for i in range(int(dim[1])):
	yvar[i]=((1-float(dim[1]))/2+i)*spacing
tvar[0]=0 
hvar[:] = S
Hvar[:] = H
Bvar[:] = B
latvar[:] = lat
lonvar[:] = lon
Accvar[:] = acc
ncfile.close()
