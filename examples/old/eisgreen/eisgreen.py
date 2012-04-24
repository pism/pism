#!/usr/bin/env python

#import Numeric
import sys
import getopt
import time
from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

GRID_FILE = 'grid20-EISMINT'
SUAQ_FILE = 'suaq20-EISMINT'
WRIT_FILE = 'eis_green20.nc'

# These values should all be outside the valid range so that generic
# applications will treat them as missing. (NUG: Attributes)
topg_fill_value = -9999.0
topg_valid_min = -5000.0

##### command line arguments #####

try:
  opts, args = getopt.getopt(sys.argv[1:], "g:p:", ["grid=","prefix="])
  for opt, arg in opts:
    if opt in ("-g", "--grid"):
      GRID_FILE = "grid" + arg + "-EISMINT"
      SUAQ_FILE = "suaq" + arg + "-EISMINT"
      WRIT_FILE = "eis_green" + arg + ".nc"
    if opt in ("-p", "--prefix"):
      GRID_FILE = arg + GRID_FILE
      SUAQ_FILE = arg + SUAQ_FILE
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

##### grid20-EISMINT #####

# open the data file and begin reading it
try:
  print "reading grid data from ",GRID_FILE
  input=open(GRID_FILE, 'r')
except IOError:
  print 'ERROR: File: ' + GRID_FILE + ' could not be found.'
  sys.exit(2)

dim=[]
for num in input.readline().split():
        dim.append(num)
input.readline() # get rid of the titles

temporary = dim[0]
dim[0] = dim[1]
dim[1] = temporary

lat = zeros((1, int(dim[0]), int(dim[1])), float32)
lon = zeros((1, int(dim[0]), int(dim[1])), float32)

latcount=0;
loncount=0
x=0
# read in all the data
for line in input.read().split():
  if x%4 == 2: # surface elevation
    lat[0, latcount/int(dim[1]), latcount%int(dim[1])]=float(line)
    latcount = latcount + 1
  elif x%4 == 3: # thickness
    lon[0, loncount/int(dim[1]), loncount%int(dim[1])]=float(line)
    loncount = loncount + 1
  x = x+1

# done reading from data file
input.close()

##### suaq20-EISMINT #####

# open the data file and begin reading it
try:
  print "reading thickness, bed elevation, and accumulation data from ",SUAQ_FILE
  input=open(SUAQ_FILE, 'r')
except IOError:
  print 'ERROR: File: ' + SUAQ_FILE + ' could not be found.'
  sys.exit(2)

dim=[]
for num in input.readline().split():
	dim.append(num)
input.readline() # get rid of the titles

temporary = dim[0]
dim[0] = dim[1]
dim[1] = temporary

S   = zeros( (int(dim[0]), int(dim[1])), float32)
H   = zeros( (int(dim[0]), int(dim[1])), float32)
B   = zeros( (int(dim[0]), int(dim[1])), float32)
acc = zeros( (int(dim[0]), int(dim[1])), float32)

ghf = 0.050 * ones( (int(dim[0]), int(dim[1])), float32)  # geothermal flux = 50 mW m-2

Scount=0;
Hcount=0
Bcount=0
Acccount=0
x=0
# read in all the data
for line in input.read().split():
	if x%4 == 0: # surface elevation
		S[Scount/int(dim[1]), Scount%int(dim[1])]=float(line)
		Scount = Scount + 1
	elif x%4 == 1: # thickness
		H[Hcount/int(dim[1]), Hcount%int(dim[1])]=float(line)
		Hcount = Hcount + 1
	elif x%4 == 2: # bedrock
		B[Bcount/int(dim[1]), Bcount%int(dim[1])]=float(line)
		Bcount = Bcount + 1
	else: # accumulation (m/a -> m/s)
		acc[Acccount/int(dim[1]), Acccount%int(dim[1])]=float(line)
		Acccount = Acccount + 1
	x = x+1

# done reading from data file
input.close()
print "Total Values Read: "+str(x)

# replace zero (used to represent missing values in the input file) with a
# value outside the valid range
putmask(B, B == 0, topg_fill_value)

# ready to write NetCDF file
ncfile = NC(WRIT_FILE, 'w',format='NETCDF3_CLASSIC')

# set global attributes
ncfile.Conventions = 'CF-1.4'
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
ncfile.history = historystr

# define the dimensions
xdim = ncfile.createDimension('x', int(dim[1]))
ydim = ncfile.createDimension('y', int(dim[0]))

# define the variables
polarVar = ncfile.createVariable('mapping', 'i4')
xvar = ncfile.createVariable('x', 'f8', dimensions=('x',))
yvar = ncfile.createVariable('y', 'f8', dimensions=('y',))
lonvar = ncfile.createVariable('lon', 'f4', dimensions=('y', 'x'))
latvar = ncfile.createVariable('lat', 'f4', dimensions=('y', 'x'))
hvar = ncfile.createVariable('usurf', 'f4', dimensions=('y', 'x'))
thkvar = ncfile.createVariable('thk', 'f4', dimensions=('y', 'x'))
bedvar = ncfile.createVariable('topg', 'f4', dimensions=('y', 'x'), fill_value=topg_fill_value) 
accvar = ncfile.createVariable('precip', 'f4', dimensions=('y', 'x'))
bheatflxvar = ncfile.createVariable('bheatflx', 'f4', dimensions=('y', 'x'))

# set the attributes of the variables
polarVar.grid_mapping_name = 'polar_stereographic'
polarVar.straight_vertical_longitude_from_pole = -41.1376
polarVar.latitude_of_projection_origin = 71.6468
polarVar.standard_parallel = 71

xvar.axis = 'X'
xvar.long_name = 'x-coordinate in Cartesian system'
xvar.standard_name = 'projection_x_coordinate'
xvar.units = 'm'

yvar.axis = 'Y'
yvar.long_name = 'y-coordinate in Cartesian system'
yvar.standard_name = 'projection_y_coordinate'
yvar.units = 'm'

lonvar.long_name = 'longitude'
lonvar.standard_name = 'longitude'
lonvar.units = 'degrees_east'

latvar.long_name = 'latitude'
latvar.standard_name = 'latitude'
latvar.units = 'degrees_north'

hvar.long_name = 'ice upper surface elevation'
hvar.standard_name = 'surface_altitude'
hvar.units = 'm'

thkvar.long_name = 'land ice thickness'
thkvar.standard_name = 'land_ice_thickness'
thkvar.units = 'm'

bedvar.long_name = 'bedrock surface elevation'
bedvar.standard_name = 'bedrock_altitude'
bedvar.units = 'm'
bedvar.valid_min = topg_valid_min

accvar.long_name = 'mean annual ice-equivalent precipitation rate'
accvar.units = 'm year-1'

bheatflxvar.long_name = 'upward geothermal flux at bedrock surface'
bheatflxvar.units = 'W m-2'

# write the data to the NetCDF file
spacing = float(dim[2])*1000
for i in range(int(dim[0])):
	yvar[i]=((1-float(dim[0]))/2+i)*spacing
for i in range(int(dim[1])):
	xvar[i]=((1-float(dim[1]))/2+i)*spacing
latvar[:] = lat
lonvar[:] = lon
hvar[:] = S
thkvar[:] = H
bedvar[:] = B
accvar[:] = acc
bheatflxvar[:] = ghf
ncfile.close()
print "NetCDF file ",WRIT_FILE," created"

