#! /usr/bin/env python

import os
import numpy
import Scientific.IO.NetCDF
import Numeric
from Scientific.IO.NetCDF import NetCDFFile
from Numeric import zeros

SECPERA = 3.1556926e7

##### grid20-EISMINT #####

# open the data file and begin reading it
input=open('grid20-EISMINT', 'r')

dim=[]
for num in input.readline().split():
        dim.append(num)
input.readline() # get rid of the titles

lat = zeros((int(dim[0]), int(dim[1])), Numeric.Float32)
lon = zeros((int(dim[0]), int(dim[1])), Numeric.Float32)

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
input=open('suaq20-EISMINT', 'r')

dim=[]
for num in input.readline().split():
	dim.append(num)
input.readline() # get rid of the titles

S = zeros((int(dim[0]), int(dim[1])), Numeric.Float32)
H = zeros((int(dim[0]), int(dim[1])), Numeric.Float32)
B = zeros((int(dim[0]), int(dim[1])), Numeric.Float32)
acc = zeros((int(dim[0]), int(dim[1])), Numeric.Float32)

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
ncfile = NetCDFFile('eis_green.nc', 'w')

# define the dimensions
ncfile.createDimension('x', int(dim[0]))
ncfile.createDimension('y', int(dim[1]))
ncfile.createDimension('t', None)

# define the variables
xvar = ncfile.createVariable('x', 'f', ('x',))
yvar = ncfile.createVariable('y', 'f', ('y',))
tvar = ncfile.createVariable('t', 'f', ('t',))
lonvar = ncfile.createVariable('lon', 'f', ('x', 'y'))
latvar = ncfile.createVariable('lat', 'f', ('x', 'y'))
hvar = ncfile.createVariable('h', 'f', ('x','y'))
Hvar = ncfile.createVariable('H', 'f', ('x','y'))
Bvar = ncfile.createVariable('b', 'f', ('x','y')) 
Accvar = ncfile.createVariable('accum', 'f', ('x','y'))

# set the attributes of the variables
setattr(xvar, 'axis', 'X')
setattr(xvar, 'long_name', 'x-coordinate in Cartesian system')
setattr(xvar, 'standard_name', 'projection_x_coordinate')
setattr(xvar, 'units', 'm')

setattr(yvar, 'axis', 'Y')
setattr(yvar, 'long_name', 'y-coordinate in Cartesian system')
setattr(yvar, 'standard_name', 'projection_y_coordinate')
setattr(yvar, 'units', 'm')

setattr(tvar, 'units', 'seconds since 1986-01-28 00:00:00 (NEED CHANGE)')

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

setattr(Accvar, 'long_name', 'mean ice equivalent accumulation rate')
setattr(Accvar, 'standard_name', 'land_ice_surface_specific_mass_balance')
setattr(Accvar, 'units', 'm y-1')

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
