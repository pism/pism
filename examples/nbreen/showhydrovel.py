#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *

import sys
import argparse

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

parser = argparse.ArgumentParser(
    description='show quiver for the subglacial water velocity (or flux) field from a PISM file')
parser.add_argument('filename',
                    help='file from which to get  V = bwatvel[2]  (and  W = bwat  for flux)')
parser.add_argument('-b', type=float, default=-1.0,
                    help='upper bound on speed; -b 100 shows all speeds above 100 as 100')
parser.add_argument('-c', type=float, default=-1.0,
                    help='arrow crop size; -c 0.1 shortens arrows longer than 0.1 * speed.max()')
# e.g. ./showhydrovel.py -c 0.001 -q extras_nbreen_y0.25_250m_routing.nc -b 20
parser.add_argument('-d', type=int, default=-1,
                    help='index of frame (default: last frame which is D=-1)')
parser.add_argument('-q', action='store_true',
                    help='show advective flux  q = V W  instead of velocity')
parser.add_argument('-s', action='store_true',
                    help='show second figure with pcolor on components of velocity (or flux)')
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
    print("ERROR: can't read from file ...")
    sys.exit(1)

print("  reading x,y axes from file %s ..." % (args.filename))
if args.t:
    xvar = nc.variables["y"]  # note x-y swap
    yvar = nc.variables["x"]
else:
    xvar = nc.variables["x"]
    yvar = nc.variables["y"]
x = asarray(squeeze(xvar[:]))
y = asarray(squeeze(yvar[:]))

print("  reading 'bwatvel[2]' field from file %s ..." % (args.filename))
try:
    velx = nc.variables["bwatvel[0]"]
except:
    print("ERROR: variable 'bwatvel[0]' not found ...")
    sys.exit(2)
try:
    vely = nc.variables["bwatvel[1]"]
except:
    print("ERROR: variable 'bwatvel[1]' not found ...")
    sys.exit(3)

if args.q:
    try:
        bwat = nc.variables["bwat"]
    except:
        print("ERROR: variable 'bwat' not found ...")
        sys.exit(6)

if args.d >= 0:
    if shape(velx)[0] <= args.d:
        print("ERROR: frame %d not available in variable velx" % args.d)
        sys.exit(3)
    if shape(vely)[0] <= args.d:
        print("ERROR: frame %d not available in variable vely" % args.d)
        sys.exit(4)
    print("  reading frame %d of %d frames" % (args.d, shape(velx)[0]))
else:
    args.d = -1
    print("  reading last frame of %d frames" % (shape(velx)[0]))

units = "m hr-1"  # FIXME: make this merely the default scale?
scale = 3.1556926e7 / 3600.0
velx = asarray(squeeze(velx[args.d, :, :])).transpose() / scale
vely = asarray(squeeze(vely[args.d, :, :])).transpose() / scale

if args.q:
    bwat = asarray(squeeze(bwat[args.d, :, :])).transpose()
    velx = velx * bwat
    vely = vely * bwat
    units = "m2 hr-1"  # FIXME: adjust units?

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
    print("  generating pcolor() image of velocity (or flux) components in figure(2) ...")
    for j in [1, 2]:
        if j == 1:
            data = velx
            name = "velx"
        else:
            data = vely
            name = "vely"
        print("  %s stats:\n    min = %9.3f %s,  max = %9.3f %s,  av = %8.3f %s" %
              (name, data.min(), units, data.max(), units, data.sum() / (x.size * y.size), units))
        subplot(1, 2, j)
        pcolor(x / 1000.0, y / 1000.0, data, vmin=data.min(), vmax=data.max())
        colorbar()
        gca().set_aspect('equal')
        gca().autoscale(tight=True)
        xlabel('x  (km)')
        ylabel('y  (km)')

speed = sqrt(velx * velx + vely * vely)

plotvelx = velx.copy()
plotvely = vely.copy()
if args.c > 0.0:
    crop = (speed > args.c * speed.max())
    plotvelx[crop] = args.c * speed.max() * velx[crop] / speed[crop]
    plotvely[crop] = args.c * speed.max() * vely[crop] / speed[crop]

if args.b > 0.0:
    speed[speed > args.b] = args.b

figure(1)
quiver(x / 1000.0, y / 1000.0, plotvelx, plotvely, speed)
colorbar()
gca().set_aspect('equal')
gca().autoscale(tight=True)
xlabel('x  (km)')
ylabel('y  (km)')

if args.q:
    print("  maximum water flux magnitude = %8.3f %s" % (speed.max(), units))
    titlestr = "water flux in %s" % units
else:
    print("  maximum water speed = %8.3f %s = %6.3f %s" %
          (speed.max(), units, speed.max() / 3600.0, 'm s-1'))  # assumes units is m hr-1
    titlestr = "water velocity in %s" % units
title(titlestr)

show()

print("  done.")
