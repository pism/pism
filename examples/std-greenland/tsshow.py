#!/usr/bin/env python

"""
Creates graph of modeled time series from multiple files.  See figures
in the "Getting Started" section of the User's Manual.
"""

# example from grid-sequencing example:
#   $ ./tsshow.py ice_volume_glacierized ice_volume_glacierized-gridseq.png ts_g20km_10ka_hy.nc '20 km for 10 ka' ts_g10km_gridseq.nc '10 km for 2 ka' ts_g5km_gridseq.nc '5 km for 200 a'
# example from paramstudy/:
#   $ ../tsshow.py ice_volume_glacierized ice_volume_glacierized-param.png ts_p10km_q0.1_e1.nc '(0.1,1)' ts_p10km_q0.5_e1.nc '(0.5,1)' ts_p10km_q1.0_e1.nc '(1.0,1)' ts_p10km_q0.1_e3.nc '(0.1,3)' ts_p10km_q0.5_e3.nc '(0.5,3)' ts_p10km_q1.0_e3.nc '(1.0,3)' ts_p10km_q0.1_e6.nc '(0.1,6)' ts_p10km_q0.5_e6.nc '(0.5,6)' ts_p10km_q1.0_e6.nc '(1.0,6)'

from numpy import *
import pylab as plt
import sys
try:
    import netCDF4 as netCDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

NC = netCDF.Dataset

if len(sys.argv) < 5:
    print("tsshow.py ERROR: at least 4 arguments needed")
    print("usage:")
    print()
    print("   $ python tsshow.py FIELD OUTIMAGE TSFILE1 LABEL1 ... TSFILEn LABELn")
    print()
    print("where strings LABEL1 ... LABELn go in the legend")
    print("example:")
    print("   $ python tsshow.py ice_volume_glacierized foo.png ts_g20km.nc '20 km' ts_g10km.nc '10 km'")
    sys.exit(1)

field = sys.argv[1]
outimage = sys.argv[2]

legloc = 'lower right'

secpera = 31556926.0
vfactor = 1.0e6 * 1.0e9

n = (len(sys.argv) - 3) / 2
labels = []
plt.figure(figsize=(9, 4))

style = ['b-',  'g-',  'r-',  'c-',  'm-',  'y-',  'k-',
         'b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for k in range(n):

    tsfile = sys.argv[2 * k + 3]
    labels.append(sys.argv[2 * k + 4])

    try:
        ncfile = NC(tsfile, "r")
    except:
        print("ERROR: can't read from file %s ..." % tsfile)
        sys.exit(2)
    t = ncfile.variables["time"][:] / secpera
    var = ncfile.variables[field][:]
    ncfile.close()
    print("read variable '%s' from time-series file '%s' ..." % (field, tsfile))

    plt.plot(t, var / vfactor, style[k], linewidth=2.5)  # automatic colors; default order
    # blue, green, red, cyan, magenta, ...
    plt.hold(True)

plt.hold(False)
plt.legend(labels, loc=legloc)
plt.xlabel("t (years)", size=16)
plt.ylabel("%s ($10^6$ km$^3$)" % field, size=16)
plt.grid(True)
print("saving image to file '%s' ..." % outimage)
# plt.show()
plt.savefig(outimage, bbox_inches='tight')
