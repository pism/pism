#!/usr/bin/env python
"""
Creates graph of the modeled time series from two files.  Such figures
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

if len(sys.argv) != 7:
    print "tsboth.py ERROR: 6 arguments needed"
    print "usage:"
    print
    print "   $ python tsboth.py FIELD TSFILE1 TSFILE2 OUTIMAGE LABEL1 LABEL2"
    print
    print "where strings LABEL1 and LABEL2 go in the legend"
    print "example:"
    print "   $ python tsboth.py ivol ts_g20km.nc ts_g10km.nc foo.png '20 km' '10 km'"
    sys.exit(1)
field = sys.argv[1]
tsfile1 = sys.argv[2]
tsfile2 = sys.argv[3]
outimage = sys.argv[4]
label1 = sys.argv[5]
label2 = sys.argv[6]

secpera = 31556926.0

try:
  nc20km = NC(tsfile1, "r")
except:
  print "ERROR: can't read from file %s ..." % tsfile1
  sys.exit(2)
t20km = nc20km.variables["time"][:]
t20km = t20km / secpera
ivol20km = nc20km.variables[field][:]
nc20km.close()
print "read variable '%s' from time-series file '%s' ..." % (field,tsfile1)

try:
  nc10km = NC(tsfile2, "r")
except:
  print "ERROR: can't read from file %s ..." % tsfile2
  sys.exit(3)
t10km = nc10km.variables["time"][:]
t10km = t10km / secpera
ivol10km = nc10km.variables[field][:]
nc10km.close()
print "read variable '%s' from time-series file '%s' ..." % (field,tsfile2)

vfactor = 1.0e6 * 1.0e9

plt.figure(figsize=(9,4))
plt.plot(t20km, ivol20km / vfactor, 'b', linewidth=2.0)
plt.hold(True)
plt.plot(t10km, ivol10km / vfactor, 'r', linewidth=2.5)
plt.hold(False)
plt.legend((label1,label2), loc='lower right')
plt.xlabel("t (years)", size=16)
plt.ylabel("%s ($10^6$ km$^3$)" % field, size=16)
plt.grid(True)

print "saving image to file '%s' ..." % outimage
plt.savefig(outimage, bbox_inches='tight')

