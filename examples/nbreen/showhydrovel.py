#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *

import sys
import argparse

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC


parser = argparse.ArgumentParser( \
    description='show quiver for the subglacial water velocity field (bwatvel) from a PISM file')
parser.add_argument('filename',
                    help='file from which to get bwatvel')
parser.add_argument('-d', type=int, default=-1,
                    help='index of frame (defaults to last frame D=-1)')
parser.add_argument('-s', action='store_true',
                    help='show second figure with pcolor on components of velocity')
parser.add_argument('-t', action='store_true',
                    help='transpose the x,y axes')
parser.add_argument('-x', action='store_true',
                    help='reverse order on x-axis')
parser.add_argument('-y', action='store_true',
                    help='reverse order on y-axis')

args = parser.parse_args()

try:
    nc = NC(args.filename, 'r')
except:
    print "ERROR: can't read from file ..."
    sys.exit(1)

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
    print "ERROR: variable 'bwatvel[0]' not found ..."
    sys.exit(2)
try:
    vely = nc.variables["bwatvel[1]"]
except:
    print "ERROR: variable 'bwatvel[1]' not found ..."
    sys.exit(3)

if args.d >= 0:
   if shape(velx)[0] <= args.d:
       print "ERROR: frame %d not available in variable velx" % args.d
       sys.exit(3)
   if shape(vely)[0] <= args.d:
       print "ERROR: frame %d not available in variable vely" % args.d
       sys.exit(4)
   print "  reading frame %d from velocity with %d frames" % (args.d, shape(velx)[0])
else:
   args.d = -1
   print "  reading last frame from velocity with %d frames" % (shape(velx)[0])

units = "m hr-1"  #FIXME: make this merely the default scale?
scale = 3.1556926e7 / 3600.0
velx = asarray(squeeze(velx[args.d,:,:])).transpose() / scale
vely = asarray(squeeze(vely[args.d,:,:])).transpose() / scale
nc.close()

if args.t:
    xytmp = velx.copy()
    velx = vely.transpose()
    vely = xytmp.transpose()
if args.x:
    x = x[::-1]
    velx = -velx
if args.y:
    y = y[::-1]
    vely = -vely

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
        print "  %s stats:\n    min = %9.3f %s,  max = %9.3f %s,  av = %8.3f %s" % \
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

print "  done."
