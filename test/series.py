#! /usr/bin/env python
## This script takes the standard out from pismr, pisms, pgrn, or pismv.  It extracts
## time series for the quantities in the summary, and saves these in a NetCDF file.
## For example,
##   $ pisms -eisII A -Mx 61 -My 61 -Mz 201 -y 5000 >> eisIIA.out
##   $ series.py -f eisIIA.out -o eisIIA_series.nc
##   $ ncview eisIIA_series.nc
## For *temperature-dependent* SIA pismv output (tests F,G,K) run as above.
## For *isothermal* SIA pismv output (tests A,B,C,D,E,L) add the option -i:
##   $ pismv -test C >> testC.out
##   $ series.py -i -f testC.out -o testC_series.nc
## CURRENT LIMITATION: doesn't work for pismd output

from numpy import *
from pycdf import *
import getopt
import sys


# defaults
IN_FILE = 'foo.txt'
OUT_FILE = 'series_out.nc'

##### process command line arguments #####
infilename = IN_FILE
outfilename = OUT_FILE
ipismvout = False
try:
  opts, args = getopt.getopt(sys.argv[1:], "i:f:o:",
                             ["isopismv=", "file=", "out="])
  for opt, arg in opts:
    if opt in ("-f", "--file"):
      infilename = arg
    elif opt in ("-o", "--out"):
      outfilename = arg
    elif opt in ("-i", "--isopismv"):
      ipismvout = True
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

# read file into an array of time and of vol
print "reading PISM (pismr, pisms, pgrn, or pismv) standard output from"
print "    ",infilename
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
  pos1 = myline.find(' (+ ')  # first string which marks a summary line
  pos2 = myline.find(']):')   # second string which marks a summary line
  posYEAR = myline.find('YEAR')
  if ((pos1 >= 13) and (pos2 >= 25) and (posYEAR < 0)):  # clearly is a summary line
    tokens = myline.split()
    year.append(float(tokens[1]))
    dtstr = tokens[3]
    bracketpos = dtstr.find('[')
    step.append(float(dtstr[:bracketpos]))
    off = 0
    if (tokens[4].find(']):') >= 0):
      off = 1
    vol.append(float(tokens[4 + off]))
    area.append(float(tokens[5 + off]))
    if (ipismvout):
      thick0.append(float(tokens[6 + off]))
    else:
      meltf.append(float(tokens[6 + off]))
      thick0.append(float(tokens[7 + off]))
      temp0.append(float(tokens[8 + off]))
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
thick0var = ncfile.def_var('thick0', NC.FLOAT, (timedim,))
setattr(thick0var, 'units', 'm')
setattr(thick0var, 'interpolation', 'linear')
if (ipismvout == False):
  meltfvar = ncfile.def_var('meltf', NC.FLOAT, (timedim,))
  setattr(meltfvar, 'interpolation', 'linear')
  temp0var = ncfile.def_var('temp0', NC.FLOAT, (timedim,))
  setattr(temp0var, 'units', 'K')
  setattr(temp0var, 'interpolation', 'linear')
# write data 
yearvar[:] = year
stepvar[:] = step
volvar[:] = vol
areavar[:] = area
thick0var[:] = thick0
if (ipismvout == False):
  meltfvar[:] = meltf
  temp0var[:] = temp0
# close
ncfile.close()
print "time series for vol, area, etc. written into NetCDF file"
print "    ", outfilename

