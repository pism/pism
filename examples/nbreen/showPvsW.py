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

parser = argparse.ArgumentParser(description='show graph of P(W) from a PISM run')
parser.add_argument('filename',
                    help='file from which to get  P = bwprel  and  W = bwat')
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

if args.d >= 0:
   if np.shape(bwprel)[0] <= args.d:
       print "ERROR: frame %d not available in variable bwprel" % args.d
       sys.exit(4)
   if np.shape(bwat)[0] <= args.d:
       print "ERROR: frame %d not available in variable bwat" % args.d
       sys.exit(5)
   print "  reading frame %d of %d frames" % (args.d, shape(bwat)[0])
else:
   args.d = -1
   print "  reading last frame of %d frames" % (np.shape(bwat)[0])

# get last frame
bwat   = np.asarray(np.squeeze(bwat[args.d,:,:])).flatten()
bwprel = np.asarray(np.squeeze(bwprel[args.d,:,:])).flatten()

bwprel[bwprel>1.0] = 1.0
bwprel[bwprel<0.0] = 0.0

nc.close()

figure(1)
plot(bwat,bwprel,'ko')
#gca().set_aspect('equal')
#gca().autoscale(tight=True)
xlabel('water thickness W  (m)')
ylabel('pressure as a fraction of overburden')

show()

print "  done."
