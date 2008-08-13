#! /usr/bin/env python

from numpy import loadtxt, concatenate
from pycdf import CDF, NC
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
ncfile = CDF(DSL_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()

# set global attributes
historysep = ' '
historystr = asctime() + ': ' + historysep.join(argv) + '\n'
setattr(ncfile, 'history', historystr)

# define time dimension, then time variable
Stdim = ncfile.def_dim('t', NC.UNLIMITED)
Stvar = ncfile.def_var('t', NC.FLOAT, (Stdim,))
setattr(Stvar, 'units', 'years before present')

# define climate data variables and attributes
d18Oseavar = ncfile.def_var('delta_18_O', NC.FLOAT, (Stdim,))
setattr(d18Oseavar, 'units', 'normalized O-18') # see specmap_readme.txt
setattr(d18Oseavar, 'long_name', 'change in oxygen isotope ratio (18^O to 16^O) from present')
setattr(d18Oseavar, 'interpolation', 'linear')
dSeavar = ncfile.def_var('delta_sea_level', NC.FLOAT, (Stdim,))
setattr(dSeavar, 'units', 'm')
setattr(dSeavar, 'long_name', 'change in sea level from present')
setattr(dSeavar, 'interpolation', 'linear')

# write data (reversing order)
DATALEN=len(years_sea)
Stvar[:DATALEN] = years_sea[::-1].tolist()
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
# compute temperature change from \delta 18^O; see C. Ritz, 1997. "EISMINT 
#   intercomparison experiment comparison of existing Greenland models"
dT = 1.5 * (d18O + 35.27)

# open the nc file to write to
ncfile = CDF(DT_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()

# set global attribute
historysep = ' '
historystr = asctime() + ': ' + historysep.join(argv) + '\n'
setattr(ncfile, 'history', historystr)

# define time dimension, variable
Ttdim = ncfile.def_dim('t', NC.UNLIMITED)
Ttvar = ncfile.def_var('t', NC.FLOAT, (Ttdim,))
setattr(Ttvar, 'units', 'years before present')

# define climate data variables and attributes
d18Ovar = ncfile.def_var('delta_18_O', NC.FLOAT, (Ttdim,))
setattr(d18Ovar, 'units', 'per mil relative to the SMOW standard') # see grip18o.readme
setattr(d18Ovar, 'long_name', 'change in isotope ratio 18^O to 16^O (oxygen) from present')
setattr(d18Ovar, 'interpolation', 'constant_piecewise_forward')
dTvar = ncfile.def_var('delta_T', NC.FLOAT, (Ttdim,))
setattr(dTvar, 'units', 'degrees C')
setattr(dTvar, 'long_name', 'change in surface temperature from present')
setattr(dTvar, 'interpolation', 'constant_piecewise_forward')

# data is redundant in pairs, except for first line which needs duplication
y_bp = concatenate(([y_bp[0]],y_bp))
y_bp = y_bp[::2]
d18O = concatenate(([d18O[0]],d18O))
d18O = d18O[::2]
dT = concatenate(([dT[0]],dT))
dT = dT[::2]

# write data
DATALEN=len(y_bp)
Ttvar[:DATALEN] = y_bp.tolist()
d18Ovar[:DATALEN] = d18O.tolist()
dTvar[:DATALEN] = dT.tolist()

# finish
ncfile.close()
print "NetCDF file ",DT_FILE," created with delta_18_O and delta_T time series"

