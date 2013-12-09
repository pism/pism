#!/usr/bin/env python

"""
Creates graph of modeled time series from multiple files.  See figures
in the "Getting Started" section of the User's Manual.
"""

from numpy import *
import pylab as plt
import sys
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
NC = netCDF.Dataset

if len(sys.argv) < 5:
    print "tsshow.py ERROR: at least 4 arguments needed"
    print "usage:"
    print
    print "   $ python tsshow.py FIELD OUTIMAGE TSFILE1 LABEL1 ... TSFILEn LABELn"
    print
    print "where strings LABEL1 ... LABELn go in the legend"
    print "example:"
    print "   $ python tsshow.py ivol foo.png ts_g20km.nc '20 km' ts_g10km.nc '10 km'"
    sys.exit(1)

field = sys.argv[1]
outimage = sys.argv[2]

legloc='lower right'

secpera = 31556926.0
vfactor = 1.0e6 * 1.0e9

n = (len(sys.argv) - 3) / 2
labels = []
plt.figure(figsize=(9,4))

style = ['b-',  'g-',  'r-',  'c-',  'm-',  'y-',  'k-',
         'b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for k in range(n):

  tsfile = sys.argv[2*k+3]
  labels.append(sys.argv[2*k+4])

  try:
    ncfile = NC(tsfile, "r")
  except:
    print "ERROR: can't read from file %s ..." % tsfile
    sys.exit(2)
  t = ncfile.variables["time"][:] / secpera
  var = ncfile.variables[field][:]
  ncfile.close()
  print "read variable '%s' from time-series file '%s' ..." % (field,tsfile)

  plt.plot(t, var / vfactor, style[k], linewidth=2.5)  # automatic colors; default order
                                             # blue, green, red, cyan, magenta, ...
  plt.hold(True)

plt.hold(False)
plt.legend(labels, loc=legloc)
plt.xlabel("t (years)", size=16)
plt.ylabel("%s ($10^6$ km$^3$)" % field, size=16)
plt.grid(True)
print "saving image to file '%s' ..." % outimage
#plt.show()
plt.savefig(outimage, bbox_inches='tight')
