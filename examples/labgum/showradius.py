#!/usr/bin/env python

# Copyright (C) 2013 the PISM Authors

# This script creates a graph of the modeled margin radius time series
# by using the iarea variable.  Compare Figure 4(a) in Sayag & Worster (2012).
# Try something like:
#    ./showradius.py -o foo.png ts_lab51.nc

from numpy import *
import pylab as plt
from sys import exit
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF

#   $ ./showradius.py -o foo.png ts_lab51.nc

from optparse import OptionParser

parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "A script to show margin radius time series a NetCDF file with 'time' dimension and 'time' and 'iarea' variables."
parser.add_option("-o", "--output_file", dest="outfile",
                  help="output file name",default='foo.png')

(options, args) = parser.parse_args()

nts = len(args)
if (nts<1) | (nts>1):
  print "needs exactly one input file ... EXITING"
  exit(-1)

nc = netCDF.Dataset(args[0], "r")

t = nc.variables["time"][:]
iarea = nc.variables["iarea"][:]
nc.close()

plt.figure(figsize=(12,6))
plt.loglog(t[t>2], sqrt(iarea[t>2]/pi) * 100.0, 'b', linewidth=2.0)  # after t=2s, and in cm
#plt.legend(('Q_0','Regression'), loc='upper left')
plt.xticks([1.0, 10.0, 100.0, 1000.0])
plt.yticks([1.0, 10.0])
plt.axis([1.0, 1000.0, 1.0, 25.0])
plt.xlabel("t  (s)", size=14)
plt.ylabel(r"$r_N$  (cm)", size=14)
plt.grid(True)

#plt.show()

plt.savefig(options.outfile)

