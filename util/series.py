#! /usr/bin/env python
## This script takes the standard out from a PISM evolution run executable, e.g.
## pismr, pisms, pgrn, or pismv.  It extracts time series for the quantities 
## in the summary, and saves these in a NetCDF file.
## For example,
##   $ pisms -eisII A -Mx 61 -My 61 -Mz 201 -y 5000 >> eisIIA.out
##   $ series.py -f eisIIA.out -o eisIIA_series.nc
##   $ ncview eisIIA_series.nc
## It works for PISM diagnostic output, but there is no useful time series
## to extract.

from numpy import *
from netCDF3 import Dataset as NC
import getopt
import sys
import time

# defaults
IN_FILE = 'foo.txt'
OUT_FILE = 'series_out.nc'

##### process command line arguments #####
infilename = IN_FILE
outfilename = OUT_FILE
ipismvout = False
try:
  opts, args = getopt.getopt(sys.argv[1:], "f:o:",
                             ["file=", "out="])
  for opt, arg in opts:
    if opt in ("-f", "--file"):
      infilename = arg
    elif opt in ("-o", "--out"):
      outfilename = arg
except getopt.GetoptError:
  print 'ERROR: incorrect command line arguments'
  sys.exit(2)

# read file into an array of time and of vol
print " reading PISM evolution run standard output from ",infilename
try:
  infile = open(infilename, 'r');
except IOError:
  print 'ERROR: file ' + infilename + ' could not be found.'
  sys.exit(2)

names=[]
units=[]
Nnames = 0
year=[]
step=[0]
vals=[]
count = 0
prototypeFound = False
unitsFound = False
ptokens = []
while True:
  myline = infile.readline()
  if not myline:
    break
  if ((myline[0] == 'P') and (myline[1] == ' ')):
    # clearly a prototype line
    tokens = myline.split()
    # expect: tokens[0] == 'P', tokens[1] == (the year): 
    #print tokens
    if (prototypeFound == True):
      # compare names and stop if different
      for j in range(Nnames):
        if names[j] != tokens[2 + j]:
          print 'ERROR: another prototype line found with conflicting name'
          print '  ("' + str(tokens[2+j]) + '" versus "' + str(names[j]) + '")'
          print 'stopping'
          sys.exit(2)
      print ' another prototype line found with identical names; continuing'
    else:
      prototypeFound = True
      Nnames = len(tokens) - 2;
      for j in range(Nnames):
        names.append(tokens[2 + j])
        vals.append([])
      print ' prototype line found; names = ' + str(names)
      #print names
  if ((myline[0] == 'U') and (myline[1] == ' ')):
    if len(units) == 0: # only read units if it is first instance
      # a units line
      tokens = myline.split()
      # expect: tokens[0] == 'U', tokens[1] == 'years'
      #print tokens
      if (Nnames == 0):
        print 'ERROR: units line found but no names defined; prototype must precede units.'
        sys.exit(2)
      for j in range(Nnames):
        units.append(tokens[2 + j])
      print ' units line found; units = ' + str(units)
  if ((myline[0] == 'S') and (myline[1] == ' ')):  # a summary line
    if (len(vals) == 0):
      print """ERROR: summary line found but no value sequences defined;\n
               prototype must precede summaries."""
      sys.exit(2)
    if (len(vals) != Nnames):
      print 'ERROR: summary line found but value sequences wrong length\n'
      sys.exit(2)
    tokens = myline.split()
    # print tokens
    yeartext = tokens[1]
    year.append(float(yeartext[0:-1]))  # remove trailing colon
    if (len(year) > 1):
      step.append(year[-1] - year[-2])
    for j in range(Nnames):
      samepos = tokens[2 + j].find('<same>')
      if samepos >= 0:
        vals[j].append(vals[j][-1])
      else:
        vals[j].append(float(tokens[2 + j]))
    count = count + 1
print ' ' + str(count) + ' summary lines read'
infile.close()

# open a NetCDF file to write to
ncfile = NC(outfilename, "w")

# add history global attribute
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(ncfile, 'history', historystr)

# always have time dimension t and step var delta_t
# define time dimension, then time variable, then attributes
NYEAR = size(year)
timedim = ncfile.createDimension('t', None)
yearvar = ncfile.createVariable('t', 'f4', dimensions=('t',))
setattr(yearvar, 'units', 'years from start of run')
stepvar = ncfile.createVariable('delta_t', 'f4', dimensions=('t',))
setattr(stepvar, 'units', 'years')
setattr(stepvar, 'interpolation', 'constant')

# define the rest of the vars
var=[]
for j in range(Nnames):
  var.append(ncfile.createVariable(names[j], 'f4', dimensions=('t',)))
  setattr(var[j], 'units', units[j])
  setattr(var[j], 'interpolation', 'linear')

# write data 
yearvar[:NYEAR] = year
stepvar[:NYEAR] = step
for j in range(Nnames):
  var[j][:NYEAR] = vals[j]

# close
ncfile.close()
print " time series written into NetCDF file ",outfilename

