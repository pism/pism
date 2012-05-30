#!/usr/bin/env python

from numpy import loadtxt, zeros
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from getopt import getopt, GetoptError
from sys import exit, argv
from time import asctime

GRIP_FILE = 'sum89-92-ss09-50yr.stp'
SPEC_FILE = 'specmap.017'
DT_FILE = 'grip_dT.nc'
DSL_FILE = 'specmap_dSL.nc'

##### command line arguments #####
try:
  opts, args = getopt(argv[1:], "p:", ["prefix="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      GRIP_FILE = arg + GRIP_FILE
      SPEC_FILE = arg + SPEC_FILE
except GetoptError:
  print 'Incorrect command line arguments'
  exit(2)


# read specmap.017
print "reading data from ",SPEC_FILE
try:
  years_sea,d18Osea = loadtxt(SPEC_FILE,skiprows=2,unpack=True)
except IOError:
  print 'ERROR: File: ' + SPEC_FILE + ' could not be found.'
  exit(2)
years_sea *= 1000.0
# compute sea level from \delta 18^O; see C. Ritz, 1997. "EISMINT 
#   intercomparison experiment comparison of existing Greenland models"
dSea = -34.83 * (d18Osea + 1.93)

# open the nc for delta Sea Level file to write to
ncfile = NC(DSL_FILE, 'w',format='NETCDF3_CLASSIC')

# set global attributes
historysep = ' '
historystr = asctime() + ': ' + historysep.join(argv) + '\n'
setattr(ncfile, 'history', historystr)

# define time dimension, then time variable
Stdim = ncfile.createDimension('t', None)
Stvar = ncfile.createVariable('t', 'f4', ('t',))
setattr(Stvar, 'units', 'years since 1-1-1')

# define climate data variables and attributes
d18Oseavar = ncfile.createVariable('delta_18_O', 'f4', ('t',))
setattr(d18Oseavar, 'units', 'normalized O-18') # see specmap_readme.txt
setattr(d18Oseavar, 'long_name', 'change in oxygen isotope ratio (18^O to 16^O) from present')
setattr(d18Oseavar, 'interpolation', 'linear')
dSeavar = ncfile.createVariable('delta_sea_level', 'f4', ('t',))
setattr(dSeavar, 'units', 'm')
setattr(dSeavar, 'long_name', 'change in sea level from present')

# write data (reversing order)
DATALEN=len(years_sea)
Stvar[:DATALEN] = ((-1) * years_sea[::-1]).tolist()
d18Oseavar[:DATALEN] = d18Osea[::-1].tolist()
dSeavar[:DATALEN] = dSea[::-1].tolist()

# finish sea level file
ncfile.close()
print "NetCDF file ",DSL_FILE," created with delta_18_O and delta_sea_level time series"


# read sum89-92-ss09-50yr.stp with GRIP data
print "reading data from ",GRIP_FILE
try:
  y_bp,d18O = loadtxt(GRIP_FILE,skiprows=9,unpack=True)
except IOError:
  print 'ERROR: File: ' + GRIP_FILE + ' could not be found.'
  exit(2)

# open the nc file to write to
ncfile = NC(DT_FILE, 'w')

# set global attribute
historysep = ' '
historystr = asctime() + ': ' + historysep.join(argv) + '\n'
setattr(ncfile, 'history', historystr)

# define time dimension, variable
Ttdim = ncfile.createDimension('t', None)
Ttvar = ncfile.createVariable('t', 'f4', ('t',))
setattr(Ttvar, 'units', 'years since 1-1-1')

# define climate data variables and attributes
d18Ovar = ncfile.createVariable('delta_18_O', 'f4', ('t',))
setattr(d18Ovar, 'units', 'per mil relative to the SMOW standard') # see grip18o.readme
setattr(d18Ovar, 'long_name', 'change in isotope ratio 18^O to 16^O (oxygen) from present')
setattr(d18Ovar, 'interpolation', 'constant_piecewise_forward')
dTvar = ncfile.createVariable('delta_T', 'f4', ('t',))
setattr(dTvar, 'units', 'Celsius')
setattr(dTvar, 'long_name', 'change in surface temperature from present')
setattr(dTvar, 'interpolation', 'constant_piecewise_forward')

# The GRIP_FILE contains 9997 records, corresponding to 4998 intervals (plus
# one extra record, which I (CK) ignored). This code computes mid-points of all
# the intervals, effectively getting rid of redundancies.
N = (len(y_bp) - 1)/2
t = zeros(N)
dOxygen = zeros(N)
for j in range(N):
  t[j] = (y_bp[2*j] + y_bp[2*j + 1]) / 2.0
  dOxygen[j] = d18O[2*j]

# compute temperature change from \delta 18^O; see C. Ritz, 1997. "EISMINT 
#   intercomparison experiment comparison of existing Greenland models"
dT = 1.5 * (dOxygen + 35.27)

# write data (reversing order)
Ttvar[:N] = ((-1) * t[::-1]).tolist()
d18Ovar[:N] = dOxygen[::-1].tolist()
dTvar[:N] = dT[::-1].tolist()

# finish
ncfile.close()
print "NetCDF file ",DT_FILE," created with delta_18_O and delta_T time series"

