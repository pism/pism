#!/usr/bin/env python

import numpy as np
from matplotlib.pyplot import *

import sys
import argparse

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 is not installed!"
    sys.exit(1)

parser = argparse.ArgumentParser(description='show scatter plot P versus W from a PISM run')
parser.add_argument('filename',
                    help='file from which to get  P = bwprel  and  W = bwat')
parser.add_argument('-c', default=None,
                    help='name of variable to use to color points in scatter plot')
parser.add_argument('-cmin', type=float, default=None,
                    help='crop color values below at this value')
parser.add_argument('-cmax', type=float, default=None,
                    help='crop color values above at this value')
parser.add_argument('-s', default=None,
                    help='name of variable to use to select whether points appear in scatter plot')
parser.add_argument('-smin', type=float, default=None,
                    help='select minimum: if -c is used then below this value (of -s var) the points will not be plotted')
parser.add_argument('-smax', type=float, default=None,
                    help='select minimum: if -c is used then above this value (of -s var) the points will not be plotted')
parser.add_argument('-d', type=int, default=-1,
                    help='index of frame (default: last frame which is D=-1)')

args = parser.parse_args()

try:
    nc = NC(args.filename, 'r')
except:
    print "ERROR: can't read from file ..."
    sys.exit(1)

print "  reading 'bwprel' field from file %s ..." % (args.filename)
try:
    bwprel = nc.variables["bwprel"]
except:
    print "ERROR: variable 'bwprel' not found ..."
    sys.exit(2)

print "  reading 'bwat' field from file %s ..." % (args.filename)
try:
    bwat = nc.variables["bwat"]
except:
    print "ERROR: variable 'bwat' not found ..."
    sys.exit(3)

if args.s != None:
  print "  reading '%s' field, for point selection, from file %s ..." % (args.s, args.filename)
  try:
      ss = nc.variables[args.s]
  except:
      print "ERROR: variable '%s' not found ..." % (args.s)
      sys.exit(5)
  if args.smin != None:
    print "    minimum value of '%s' for selection is %f" % (args.s, args.smin)
  if args.smax != None:
    print "    maximum value of '%s' for selection is %f" % (args.s, args.smax)

if args.c != None:
  print "  reading '%s' field, for point color, from file %s ..." % (args.c, args.filename)
  try:
      cc = nc.variables[args.c]
  except:
      print "ERROR: variable '%s' not found ..." % (args.c)
      sys.exit(4)
  if args.cmin != None:
    print '    color minimum value is %f' % args.cmin
  if args.cmax != None:
    print '    color maximum value is %f' % args.cmax

if args.d >= 0:
   if np.shape(bwprel)[0] <= args.d:
       print "ERROR: frame %d not available in variable bwprel" % args.d
       sys.exit(4)
   if np.shape(bwat)[0] <= args.d:
       print "ERROR: frame %d not available in variable bwat" % args.d
       sys.exit(5)
   print "  using frame %d of %d frames" % (args.d, np.shape(bwat)[0])
else:
   args.d = -1
   print "  reading last frame of %d frames" % (np.shape(bwat)[0])

bwat   = np.asarray(np.squeeze(bwat[args.d,:,:])).flatten()
bwprel = np.asarray(np.squeeze(bwprel[args.d,:,:])).flatten()

bwprel[bwprel>1.0] = 1.0
bwprel[bwprel<0.0] = 0.0

if args.c != None:
  ccc = np.asarray(np.squeeze(cc[args.d,:,:])).flatten()
  if args.cmin != None:
    ccc[ccc<args.cmin] = args.cmin
  if args.cmax != None:
    ccc[ccc>args.cmax] = args.cmax

if args.s != None:
  sss = np.asarray(np.squeeze(ss[args.d,:,:])).flatten()
  if args.smin != None:
    bwat = bwat[sss>=args.smin]
    bwprel = bwprel[sss>=args.smin]
    if args.c != None:
      ccc = ccc[sss>=args.smin]
    sss = sss[sss>=args.smin]
  if args.smax != None:
    bwat = bwat[sss<=args.smax]
    bwprel = bwprel[sss<=args.smax]
    if args.c != None:
      ccc = ccc[sss<=args.smax]
    sss = sss[sss<=args.smax]

nc.close()

figure(1)
if args.c != None:
  scatter(bwat,bwprel,c=ccc)
  colorbar()
else:
  scatter(bwat,bwprel,c='k')
#gca().set_aspect('equal')
gca().autoscale(tight=True)
xlabel('water thickness W  (m)')
ylabel('pressure as a fraction of overburden')

show()

print "  done."
