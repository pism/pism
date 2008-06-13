#! /usr/bin/env python

from numpy import *
from pycdf import *
import getopt
import sys

GRIP_FILE = 'sum89-92-ss09-50yr.stp'
SPEC_FILE = 'specmap.017'
DT_FILE = 'grip_dT.nc'
DSL_FILE = 'specmap_dSL.nc'
SPEC_LENGTH = 782

##### command line arguments #####
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:", ["prefix="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      GRIP_FILE = arg + GRIP_FILE
      SPEC_FILE = arg + SPEC_FILE
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

##### sea level from SPECMAP first  #####

# read specmap.017
print "reading data from ",SPEC_FILE
try:
  input = open(SPEC_FILE, 'r');
except IOError:
  print 'ERROR: File: ' + SPEC_FILE + ' could not be found.'
  sys.exit(2)

years_sea=[]
d18Osea=[]
dSea=[]
# ignor some lines
input.readline()
input.readline()
# take care of the fact that the data is offset
d18Osea.append(0.0)
dSea.append(0.0)
count=0;
for num in input.read().split():
  if count%2 == 0:
    years_sea.append(float(num) * 1000)
  else:
    d18Osea.append(float(num))
    # compute delta sea level from formula in EISMINT-Greenland description
    dSea.append(-34.83 * (float(num) + 1.93))
  count = count + 1
print str(count) + ' numbers read.'
d18Osea = d18Osea[:(size(d18Osea)-1)]
dSea = dSea[:(size(dSea)-1)]
years_sea.reverse()
d18Osea.reverse()
dSea.reverse()
input.close()

# open the nc for delta Sea Level file to write to
ncfile = CDF(DSL_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()
# define time dimension, then time variable, then attributes
Stdim = ncfile.def_dim('t', size(dSea))
Stvar = ncfile.def_var('t', NC.FLOAT, (Stdim,))
setattr(Stvar, 'units', 'years before present')
d18Oseavar = ncfile.def_var('delta_18_O', NC.FLOAT, (Stdim,))
setattr(d18Oseavar, 'units', 'normalized O-18') # see specmap_readme.txt
setattr(d18Oseavar, 'long_name', 'change in oxygen isotope ratio (18^O to 16^O) from present')
setattr(d18Oseavar, 'interpolation', 'linear')
dSeavar = ncfile.def_var('delta_sea_level', NC.FLOAT, (Stdim,))
setattr(dSeavar, 'units', 'm')
setattr(dSeavar, 'long_name', 'change in sea level from present')
setattr(dSeavar, 'interpolation', 'linear')
# write data and close
Stvar[:] = years_sea
d18Oseavar[:] = d18Osea
dSeavar[:] = dSea
ncfile.close()
print "NetCDF file ",DSL_FILE," created"


##### delta T (and delta O18) from GRIP next  #####
# read sum89-92-ss09-50yr.stp
dim=[]
try:
  print "reading data from ",GRIP_FILE
  input = open(GRIP_FILE, 'r')
except IOError:
  print 'ERROR: File: ' + GRIP_FILE + ' could not be found.'
  sys.exit(2)

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
    y_bp[y_count] = float(num)
  if count%4 == 1:
    d18O[y_count] = float(num)
    # compute delta T from formula in EISMINT-Greenland description
    dT[y_count] = 1.5*(d18O[y_count]+35.27)
    y_count = y_count-1
  count = count + 1
print str(count) + ' numbers read.'

# open the nc file to write to
ncfile = CDF(DT_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()
# define time dimension, then time variable, then attributes
Ttdim = ncfile.def_dim('t', int(dim[0])/2)
Ttvar = ncfile.def_var('t', NC.FLOAT, (Ttdim,))
setattr(Ttvar, 'units', 'years before present')
d18Ovar = ncfile.def_var('delta_18_O', NC.FLOAT, (Ttdim,))
setattr(d18Ovar, 'units', 'per mil relative to the SMOW standard') # see grip18o.readme
setattr(d18Ovar, 'long_name', 'change in isotope ratio 18^O to 16^O (oxygen) from present')
setattr(d18Ovar, 'interpolation', 'constant_piecewise_forward')
dTvar = ncfile.def_var('delta_T', NC.FLOAT, (Ttdim,))
setattr(dTvar, 'units', 'degrees C')
setattr(dTvar, 'long_name', 'change in surface temperature from present')
setattr(dTvar, 'interpolation', 'constant_piecewise_forward')
# write data and close
Ttvar[:] = y_bp
d18Ovar[:] = d18O
dTvar[:] = dT
ncfile.close()
print "NetCDF file ",DT_FILE," created"

