#! /usr/bin/env python

from numpy import *
from pycdf import *
import getopt
import sys

SUM8_FILE = 'sum89-92-ss09-50yr.stp'
SPEC_FILE = 'specmap.017'
CORE_FILE = 'eis_core.nc'
INTERPOLATION = 'linear'

SPEC_LENGTH = 782

##### command line arguments #####

try:
  opts, args = getopt.getopt(sys.argv[1:], "i:", ["interpolation"])
  for opt, arg in opts:
    if opt in ("-i", "--interpolation"):
      INTERPOLATION = arg
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

##### specmap.017 #####

input = open(SPEC_FILE, 'r');

years_sea=[]
dSea=[]
#years_sea = zeros(SPEC_LENGTH, float32)
#dSea = zeros(SPEC_LENGTH, float32)

input.readline()
input.readline()

# take care of the fact that the data is offset
dSea.append(0.0)

count=0;
for num in input.read().split():
  if count%2 == 0:
    years_sea.append(-float(num) * 1000)
  else:
    dSea.append(-34.83 * (float(num) + 1.93))
  count = count + 1

dSea = dSea[:(size(dSea)-1)]

years_sea.reverse()
dSea.reverse()

print "len(years_sea) = ",size(years_sea)
print "len(dSea) = ",size(dSea)

input.close()

##### sum89-92-ss09-50yr.stp #####

dim=[]

input = open(SUM8_FILE, 'r')

# remove headers 
for n in zeros(8):
  input.readline()

for num in input.readline().split():
  dim.append(num)

y_bp = zeros(int(dim[0])/2, float32)
d18O = zeros(int(dim[0])/2, float32)
dT = zeros(int(dim[0])/2, float32)
count = 0
y_count= int(dim[0])/2-1

for num in input.read().split():
  if count%4 == 0:
    y_bp[y_count] = -float(num)
  if count%4 == 1:
    d18O[y_count] = float(num)
    dT[y_count] = 1.5*(d18O[y_count]+35.27)
    y_count = y_count-1
  count = count + 1

print str(count) + ' numbers read.'

# open the nc file to write to
ncfile = CDF(CORE_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()


# define the one dimension
Ttdim = ncfile.def_dim('delta_T_t', int(dim[0])/2)
Stdim = ncfile.def_dim('delta_Sea_t', size(dSea))

# define variables
polarVar = ncfile.def_var('polar_stereographic', NC.INT)
Ttvar = ncfile.def_var('delta_T_t', NC.FLOAT, (Ttdim,))
Stvar = ncfile.def_var('delta_Sea_t', NC.FLOAT, (Stdim,))
d18Ovar = ncfile.def_var('delta_O_18', NC.FLOAT, (Ttdim,))
dTvar = ncfile.def_var('delta_T', NC.FLOAT, (Ttdim,))
dSeavar = ncfile.def_var('delta_Sea', NC.FLOAT, (Stdim,))

# define attributes
setattr(polarVar, 'grid_mapping_name', 'polar_stereographic')
setattr(polarVar, 'straight_vertical_longitude_from_pole', 0)
setattr(polarVar, 'straight_vertical_longitude_from_pole', 90)
setattr(polarVar, 'standard_parallel', -71)

setattr(Ttvar, 'units', 'years since 1989')

setattr(Stvar, 'units', 'years since 1989')

setattr(dTvar, 'units', 'm')
setattr(dTvar, 'long_name', 'change in surface temperature from 1986')
setattr(dTvar, 'interpolation', INTERPOLATION)

setattr(dSeavar, 'units', 'm')
setattr(dSeavar, 'long_name', 'change in sea level from 1986')
setattr(dSeavar, 'interpolation', INTERPOLATION)

# write data
Ttvar[:] = y_bp
Stvar[:] = years_sea
d18Ovar[:] = d18O
dTvar[:] = dT
dSeavar[:] = dSea

ncfile.close()
