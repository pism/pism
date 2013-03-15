#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *
from sys import exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

import argparse

parser = argparse.ArgumentParser( \
    description='show quiver for the subglacial water velocity field (bwatvel) from a PISM file')
parser.add_argument('filename',
                    help='file from which to get bwatvel')
parser.add_argument('-c', nargs=4, type=int, metavar='x',
                    help='crop to rectangle [xmin xmax ymin ymax]; in km')
parser.add_argument('-s', action='store_true',
                    help='show second figure with pcolor on components of velocity')
parser.add_argument('-t', action='store_true',
                    help='use transpose on x,y axes')

args = parser.parse_args()

try:
    nc = NC(args.filename, 'r')
except:
    print "  ERROR: can't read from file ..."
    exit(1)

#FIXME:  add option to crop to rectangle given in km

print "  reading x,y axes from file %s ..." % (args.filename)
if args.t:
    xvar = nc.variables["y"]  # note x-y swap
    yvar = nc.variables["x"]
else:
    xvar = nc.variables["x"]
    yvar = nc.variables["y"]
x = asarray(squeeze(xvar[:]))
y = asarray(squeeze(yvar[:]))

print "  reading 'bwatvel[2]' field from file %s ..." % (args.filename)
try:
    velx = nc.variables["bwatvel[0]"]
except:
    print "  ERROR variable 'bwatvel[0]' not found ..."
    exit(2)
try:
    vely = nc.variables["bwatvel[1]"]
except:
    print "  ERROR variable 'bwatvel[1]' not found ..."
    exit(3)

units = "m hr-1"  #FIXME: make this merely the default scale
scale = 3.1556926e7 / 3600.0
velx = asarray(squeeze(velx[-1,:,:])).transpose() / scale
vely = asarray(squeeze(vely[-1,:,:])).transpose() / scale

if args.t:
    xytmp = velx.copy()
    velx = vely.transpose()
    vely = xytmp.transpose()

print args.c

if args.s:
    figure(2)
    print "  generating pcolor() image of velocity components in figure(2) ..."
    for j in [1,2]:
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

figure(1)
speed = sqrt(velx*velx + vely*vely)
print "  maximum water speed = %8.3f %s" % (speed.max(),units)
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
