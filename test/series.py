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
  opts, args = getopt.getopt(sys.argv[1:], "f:o:",
                             ["file=", "out="])
  for opt, arg in opts:
    if opt in ("-f", "--file"):
      infilename = arg
    elif opt in ("-o", "--out"):
      outfilename = arg
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

# read file into an array of time and of vol
print "reading PISM evolution run standard output from ",infilename
try:
  infile = open(infilename, 'r');
except IOError:
  print 'ERROR: File: ' + infilename + ' could not be found.'
  sys.exit(2)

names=[]
units=[]
Nnames = 0
year=[]
step=[]
vals=[]
count = 0
prototypeFound = False
while True:
  myline = infile.readline()
  if not myline:
    break
  posbracketplus = myline.find('(+')  # second marker of a prototype or summary line
  if ((myline[0] == 'P') and (myline[1] == ' ')
      and (posbracketplus >= 0) and (prototypeFound == True)):
    print 'second prototype line found; stopping'
    break
  if ((myline[0] == 'P') and (myline[1] == ' ')
      and (posbracketplus >= 0) and (prototypeFound == False)):
    # clearly a prototype line
    tokens = myline.split()
    # expect: tokens[0] == 'S', tokens[1] == (the year), tokens[2] == '(', 
    #         tokens[3] == (the step), tokens[4] == '[N$])'
    #print tokens
    Nnames = len(tokens) - 5;
    for j in range(Nnames):
      names.append(tokens[5 + j])
      vals.append([])
    #print 'prototype line found, with names = '
    #print names
  if ((myline[0] == 'U') and (myline[1] == ' ')):
    # clearly a units line
    tokens = myline.split()
    # expect: tokens[0] == 'U', tokens[1] == 'years', tokens[2] == 'years'
    #print tokens
    if (Nnames == 0):
      print 'ERROR: units line found but no names defined; prototype must precede units.'
      sys.exit(2)
    for j in range(Nnames):
      units.append(tokens[3 + j])
    #print 'units line found, with units = '
    #print units
  if ((myline[0] == 'S') and (posbracketplus >= 0)):  # clearly is a summary line
    if (len(vals) == 0):
      print """ERROR: summary line found but no value sequences defined;\n
               prototype must precede summaries."""
      sys.exit(2)
    if (len(vals) != Nnames):
      print 'ERROR: summary line found but value sequences wrong length\n'
      sys.exit(2)
    tokens = myline.split()
    year.append(float(tokens[1]))
    step.append(float(tokens[3]))
    for j in range(Nnames):
      vals[j].append(float(tokens[5 + j]))
    count = count + 1
print str(count) + ' summary lines read'
infile.close()

# open a NetCDF file to write to
ncfile = CDF(outfilename, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()

# always have time dimension t and step var delta_t
# define time dimension, then time variable, then attributes
timedim = ncfile.def_dim('t', size(year))
yearvar = ncfile.def_var('t', NC.FLOAT, (timedim,))
setattr(yearvar, 'units', 'years from start of run')
stepvar = ncfile.def_var('delta_t', NC.FLOAT, (timedim,))
setattr(stepvar, 'units', 'years')
setattr(stepvar, 'interpolation', 'constant')

# define the rest of the vars
var=[]
for j in range(Nnames):
  var.append(ncfile.def_var(names[j], NC.FLOAT, (timedim,)))
  setattr(var[j], 'units', units[j])
  setattr(var[j], 'interpolation', 'linear')

# write data 
yearvar[:] = year
stepvar[:] = step
for j in range(Nnames):
  var[j][:] = vals[j]

# close
ncfile.close()
print "time series written into NetCDF file ",outfilename

