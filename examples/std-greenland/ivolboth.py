#!/usr/bin/env python
"""
This script creates a graph of the modeled volume time series
for both 20km and 10km grid calculations.  This figure
appears in the "Getting Started" section of the User's Manual.
"""

from numpy import *
import pylab as plt
import sys
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
NC = netCDF.Dataset

if len(sys.argv) < 4:
    print "ERROR: usage error"
    print "usage:"
    print "  $ python ivolboth.py ts_g20km.nc ts_g10km.nc both.png"
    sys.exit(1)

secpera = 31556926.0

try:
  nc20km = NC(sys.argv[1], "r")
except:
  print "ERROR: can't read from file %s ..." % sys.argv[1]
  sys.exit(2)
t20km = nc20km.variables["time"][:]
t20km = t20km / secpera
ivol20km = nc20km.variables["ivol"][:]
nc20km.close()

try:
  nc10km = NC(sys.argv[2], "r")
except:
  print "ERROR: can't read from file %s ..." % sys.argv[2]
  sys.exit(3)
t10km = nc10km.variables["time"][:]
t10km = t10km / secpera
ivol10km = nc10km.variables["ivol"][:]
nc10km.close()

vfactor = 1.0e6 * 1.0e9

plt.figure(figsize=(9,4))
plt.plot(t20km, ivol20km / vfactor, 'b', linewidth=2.0)
plt.hold(True)
plt.plot(t10km, ivol10km / vfactor, 'r', linewidth=2.5)
plt.hold(False)
plt.legend(('20 km','10 km'), loc='lower right')
plt.xlabel("t (years)", size=16)
plt.ylabel("volume ($10^6$ km$^3$)", size=16)
plt.grid(True)

plt.savefig(sys.argv[3], bbox_inches='tight')

