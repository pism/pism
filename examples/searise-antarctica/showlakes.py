#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *
from sys import exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

filename = "lakes100km.nc"

try:
    nc = NC(filename, 'r')
except:
    print "Note: %s was not modified." % output_filename
    exit(-1)

print "loading 'bwat' from %s ..." % filename

xvar = nc.variables["x"]
yvar = nc.variables["y"]
x = asarray(squeeze(xvar[:]))
y = asarray(squeeze(yvar[:]))

bwatvar = nc.variables["bwat"]
#print bwatvar
bwat = asarray(squeeze(bwatvar[:])).transpose()

pcolor(x/1000.0,y/1000.0,bwat)
#axis('equal')
gca().set_aspect('equal')
gca().autoscale(tight=True)
xlabel('x  (km)')
ylabel('y  (km)')

nc.close()

show()

