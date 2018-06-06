#!/usr/bin/env python

"""
Creates graph of three modeled hydrology-related time series from a single file.
"""

# example from grid-sequencing example:
#   $ ./hydro-tsshow.py foo.png ts_routing-decoupled.nc

from numpy import *
import pylab as plt
import sys
try:
    import netCDF4 as netCDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

NC = netCDF.Dataset

if len(sys.argv) < 3:
    print("hydro-tsshow.py ERROR: ... FIXME ... exiting")
    sys.exit(1)
outimage = sys.argv[1]
tsfile = sys.argv[2]

secpera = 31556926.0
scale = 10.0e3
scalestr = '$10^3$'
yaxismin = 1.0e-4  # in scale*kg/s
legloc = 'lower right'

labels = []
plt.figure(figsize=(9, 4))

print("opening file '%s' for reading ..." % tsfile)
try:
    ncfile = NC(tsfile, "r")
except:
    print("ERROR: can't open file %s for reading ..." % tsfile)
    sys.exit(2)
print("  reading 'time' variable ...")
t = ncfile.variables["time"][:] / secpera

n = 3
style = ['b-',  'g-',  'r-']
labels = ['ocean_loss', 'ice_free_land_loss', 'negative_thickness_gain']
for k in range(n):
    varname = 'hydro_' + labels[k]
    print("  reading '%s' variable ..." % varname)
    var = ncfile.variables[varname][:]
    plt.semilogy(t, var / scale, style[k], linewidth=2.5)
    plt.hold(True)

ncfile.close()

plt.hold(False)
yy = plt.getp(plt.gca(), 'ylim')
plt.setp(plt.gca(), 'ylim', (yaxismin, yy[1]))
plt.legend(labels, loc=legloc)
plt.xlabel("t (years)", size=16)
plt.ylabel("flux  (%s kg/s)" % scalestr, size=16)
plt.grid(True)

print("saving image to file '%s' ..." % outimage)
# plt.show()
plt.savefig(outimage, bbox_inches='tight')
