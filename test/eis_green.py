# Copyright (C) 2007 Jed Brown and Ed Bueler
#
# This file is part of Pism.
#
# Pism is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# Pism is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with Pism; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import os
import numpy
import Scientific.IO.NetCDF
import Numeric
from Scientific.IO.NetCDF import NetCDFFile
from Numeric import zeros

##### grid20-EISMINT #####

# open the data file and begin reading it
input=open('grid20-EISMINT', 'r')

dim=[]
for num in input.readline().split():
        dim.append(num)
input.readline() # get rid of the titles

lat = zeros((1, int(dim[0]), int(dim[1])), Numeric.Float32)
lon = zeros((1, int(dim[0]), int(dim[1])), Numeric.Float32)

latcount=0;
loncount=0
x=0
# read in all the data
for line in input.read().split():
        if x%4 == 2: # surface elevation
                lat[0, latcount%int(dim[0]), latcount/int(dim[0])]=float(line)
                latcount = latcount + 1
        elif x%4 == 3: # thickness
                lon[0,loncount%int(dim[0]), loncount/int(dim[0])]=float(line)
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

S = zeros((1, int(dim[0]), int(dim[1])), Numeric.Float32)
H = zeros((1, int(dim[0]), int(dim[1])), Numeric.Float32)
B = zeros((1, int(dim[0]), int(dim[1])), Numeric.Float32)
acc = zeros((1, int(dim[0]), int(dim[1])), Numeric.Float32)

Scount=0;
Hcount=0
Bcount=0
Acccount=0
x=0
# read in all the data
for line in input.read().split():
	if x%4 == 0: # surface elevation
		S[0, Scount%int(dim[0]), Scount/int(dim[0])]=float(line)
		Scount = Scount + 1
	elif x%4 == 1: # thickness
		H[0,Hcount%int(dim[0]), Hcount/int(dim[0])]=float(line)
		Hcount = Hcount + 1
	elif x%4 == 2: # bedrock
		B[0,Bcount%int(dim[0]), Bcount/int(dim[0])]=float(line)
		Bcount = Bcount + 1
	else: # accumulation (m/a -> m/s)
		acc[0,Acccount%int(dim[0]), Acccount/int(dim[0])]=float(line)/3.1556926e7
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
lonvar = ncfile.createVariable('lon', 'f', ('t', 'x', 'y'))
latvar = ncfile.createVariable('lat', 'f', ('t', 'x', 'y'))
spacevar = ncfile.createVariable('spacing', 'f', ('t',))
Svar = ncfile.createVariable('S', 'f', ('t', 'x','y'))
Hvar = ncfile.createVariable('H', 'f', ('t', 'x','y'))
Bvar = ncfile.createVariable('b', 'f', ('t', 'x','y')) 
Accvar = ncfile.createVariable('accum', 'f', ('t', 'x','y'))

# set the attributes of the variables
setattr(spacevar, 'long_name', 'grid spacing on the xy-plane')
setattr(spacevar, 'units', 'm')

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

setattr(Svar, 'long_name', 'surface_altitude')
setattr(Svar, 'standard_name', 'surface_altitude')
setattr(Svar, 'units', 'm')

setattr(Hvar, 'long_name', 'land_ice_thickness')
setattr(Hvar, 'standard_name', 'land_ice_thickness')
setattr(Hvar, 'units', 'm')

setattr(Bvar, 'long_name', 'bedrock_altitude')
setattr(Bvar, 'standard_name', 'bedrock_altitude')
setattr(Bvar, 'units', 'm')

setattr(Accvar, 'long_name', 'mean ice equivalent accumulation rate')
setattr(Accvar, 'standard_name', 'land_ice_surface_specific_mass_balance')
setattr(Accvar, 'units', 'm s-1')

# write the data to the NetCDF file
for i in range(int(dim[0])):
	xvar[i]=((1-int(dim[0]))/2+i)*int(dim[2])
for i in range(int(dim[1])):
	yvar[i]=((1-int(dim[1]))/2+i)*int(dim[2])
tvar[0]=0 
spacevar[0]=float(dim[2])*1000
Svar[:] = S
Hvar[:] = H
Bvar[:] = B
latvar[:] = lat
lonvar[:] = lon
Accvar[:] = acc
ncfile.close()
