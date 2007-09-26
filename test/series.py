#! /usr/bin/env python
## This script takes the standard out from pismr or pisms, extracts time series
## for the quantities in the summary, and saves these in a NetCDF file.
## For example,
##   $ pisms -eisII A -Mx 61 -My 61 -Mz 201 -y 5qq000 >> eisIIA.out
##   $ series.py -f eisIIA.out -o eisIIA_series.nc
##   $ ncview eisIIA_series.nc
##
## CURRENT LIMITATIONS:
##   1. doesn't work for pismv or pismd output
##   2. doesn't work for -verbose output

from numpy import *
from pycdf import *
import getopt
import sys


# defaults
IN_FILE = 'foo.txt'
OUT_FILE = 'series_out.nc'

##### process command line arguments #####
pref = ''
infilename = IN_FILE
outfilename = OUT_FILE
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:f:o:",
                             ["prefix=", "file=", "out="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      pref = arg
    elif opt in ("-f", "--file"):
      infilename = arg
    elif opt in ("-o", "--out"):
      outfilename = arg
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

# read file into an array of time and of vol
print "reading PISM (pisms or pismr) standard output from ",infilename
try:
  infile = open(infilename, 'r');
except IOError:
  print 'ERROR: File: ' + infilename + ' could not be found.'
  sys.exit(2)

year=[]
step=[]
vol=[]
area=[]
meltf=[]
thick0=[]
temp0=[]
count = 0
while True:
  myline = infile.readline()
  if not myline:
    break
  if myline[0] == '$':
    year.append(float(myline[5:16]))
    step.append(float(myline[19:28]))
    vol.append(float(myline[34:42]))
    area.append(float(myline[42:50]))
    meltf.append(float(myline[50:59]))
    thick0.append(float(myline[59:70]))
    temp0.append(float(myline[70:80]))
    count = count + 1
print str(count) + ' summary lines read'
infile.close()

# open a NetCDF file to write to
ncfile = CDF(outfilename, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()
# define time dimension, then time variable, then attributes
timedim = ncfile.def_dim('t', size(year))
yearvar = ncfile.def_var('t', NC.FLOAT, (timedim,))
setattr(yearvar, 'units', 'years from start of run')
# define the rest of the vars
stepvar = ncfile.def_var('delta_t', NC.FLOAT, (timedim,))
setattr(stepvar, 'units', 'years')
setattr(stepvar, 'interpolation', 'constant')
volvar = ncfile.def_var('vol', NC.FLOAT, (timedim,))
setattr(volvar, 'units', '10^6 km^3')
setattr(volvar, 'interpolation', 'linear')
areavar = ncfile.def_var('area', NC.FLOAT, (timedim,))
setattr(areavar, 'units', '10^6 km^2')
setattr(areavar, 'interpolation', 'linear')
meltfvar = ncfile.def_var('meltf', NC.FLOAT, (timedim,))
setattr(meltfvar, 'interpolation', 'linear')
thick0var = ncfile.def_var('thick0', NC.FLOAT, (timedim,))
setattr(thick0var, 'units', 'm')
setattr(thick0var, 'interpolation', 'linear')
temp0var = ncfile.def_var('temp0', NC.FLOAT, (timedim,))
setattr(temp0var, 'units', 'K')
setattr(temp0var, 'interpolation', 'linear')
# write data 
yearvar[:] = year
stepvar[:] = step
volvar[:] = vol
areavar[:] = area
meltfvar[:] = meltf
thick0var[:] = thick0
temp0var[:] = temp0
# close
ncfile.close()
print "time series for vol, area, etc. written into NetCDF file ",outfilename

