#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *
from sys import exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

from optparse import OptionParser

parser = OptionParser()

parser.usage = "%prog [options]"
parser.description = "Shows the hydrologic velocity field 'bwatvel' from a PISM-generated file."
parser.add_option("-f", dest="filename", type="string",
                  help="file from which to get bwatvel")

(opts, args) = parser.parse_args()

filename = opts.filename

#print "  attempting to read 'bwatvel[2]' field from file %s ..." % (filename)

try:
    nc = NC(filename, 'r')
except:
    print "  ERROR: can't read from file ..."
    exit(1)

xvar = nc.variables["y"]  # note x-y swap
yvar = nc.variables["x"]
x = asarray(squeeze(xvar[:]))
y = asarray(squeeze(yvar[:]))

print "  generating pcolor() image of velocity components ..."
try:
    velx = nc.variables["bwatvel[1]"]
except:
    print "  ERROR variable 'bwatvel[1]' not found ..."
    exit(2)
try:
    vely = nc.variables["bwatvel[0]"]
except:
    print "  ERROR variable 'bwatvel[0]' not found ..."
    exit(3)

scale = 3.1556926e7 / 3600.0
velx = asarray(squeeze(velx[-1,:,:])) / scale
vely = asarray(squeeze(vely[-1,:,:])) / scale

units = "m hr-1"

if False:
#figure(1)
#for j in [1,2]:
    if j == 1:
        data = velx
        name = "velx"
    else:
        data = vely
        name = "vely"
    print "  %s stats:  min = %9.3f %s,  max = %9.3f %s,  av = %8.3f %s" % \
          (name,data.min(),units,data.max(),units,data.sum()/(x.size*y.size),units)
    subplot(1,2,j)
    pcolor(x/1000.0,y/1000.0,data,vmin=data.min(),vmax=data.max())
    colorbar()
    gca().set_aspect('equal')
    gca().autoscale(tight=True)
    xlabel('x  (km)')
    ylabel('y  (km)')

figure(2)
speed = sqrt(velx*velx + vely*vely)
quiver(x/1000.0,y/1000.0,velx,vely,speed)
colorbar()
title("water velocity in meters/hour")
gca().set_aspect('equal')
gca().autoscale(tight=True)
xlabel('x  (km)')
ylabel('y  (km)')

show()

#savefig(pngfilename, dpi=300, bbox_inches='tight')
#close()

nc.close()
print "  done."
