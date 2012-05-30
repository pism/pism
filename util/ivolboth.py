#!/usr/bin/env python
"""This script creates a graph of the modeled volume time series
for both 20km and 10km grid calculations for SeaRISE-Greenland.  This figure
appears in the "Getting Started" section of the User's Manual.
"""

from numpy import *
import pylab as plt
from sys import exit
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
NC = netCDF.Dataset

# generate "whole" time series this way:
#   $ ncrcat ts_g20km_m5000a.nc ts_g20km_0.nc -o ts_g20km_whole.nc
#   $ ncrcat ts_g10km_m5000a.nc ts_g10km_0.nc -o ts_g10km_whole.nc
# run this way:
#   $ python ivolboth.py -o both.png ts_g20km_whole.nc ts_g10km_whole.nc

from optparse import OptionParser

parser = OptionParser()
parser.usage = "usage: %prog [options] FILE1 FILE2"
parser.description = "A script to show ivol time series from two files."
parser.add_option("-o", "--output_file", dest="outfile",
                  help="output file name",default='both.png')

(options, args) = parser.parse_args()
nts = len(args)
if (nts<2) | (nts>2):
  print "needs exactly two input files ... EXITING"
  exit(-1)

nc20km = NC(args[0], "r")

secpera = 31556926.0
t20km = nc20km.variables["time"][:]
t20km = t20km / secpera
ivol20km = nc20km.variables["ivol"][:]
nc20km.close()

nc10km = NC(args[1], "r")
t10km = nc10km.variables["time"][:]
t10km = t10km / secpera
ivol10km = nc10km.variables["ivol"][:]
nc10km.close()

vfactor = 1.0e6 * 1.0e9

plt.figure(figsize=(12,6))
plt.plot(t20km, ivol20km / vfactor, 'b', linewidth=2.0)
plt.hold(True)
plt.plot(t10km, ivol10km / vfactor, 'r', linewidth=2.5)
plt.hold(False)
plt.legend(('20 km','10 km'), loc='lower right')
plt.xlabel("t (years)", size=16)
plt.ylabel("volume ($10^6$ km$^3$)", size=16)
plt.grid(True)

#last=-10000.0
#t20km_last = t20km[t20km >= last]
#ivol20km_last = ivol20km[t20km >= last]
#t10km_last = t10km[t10km >= last]
#ivol10km_last = ivol10km[t10km >= last]

#axesinset = axes([0.55, 0.20, 0.33, 0.25],axisbg='w')
#plot(t20km_last,ivol20km_last / vfactor, 'b', linewidth=1.0), hold(True)
#plot(t10km_last,ivol10km_last / vfactor, 'r', linewidth=2.0), hold(False)
#setp(axesinset)

#setp(axesinset,xticks=arange(95.e3,101.e3,1.e3),
#     xticklabels=('95','96','97','98','99','100'))

#plt.show()

plt.savefig(options.outfile)

