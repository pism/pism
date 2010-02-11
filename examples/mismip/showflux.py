#!/usr/bin/env python

from pylab import *
from netCDF3 import Dataset as NC
import os
import sys
from getopt import getopt, GetoptError


infilename=""
outfilename="foo.png"
try:
  opts, args = getopt(sys.argv[1:], "f:o:", ['filename=','outname='])
  for opt, arg in opts:
    if opt in ("-f", "--filename"):
      infilename = arg
    if opt in ("-o", "--outname"):
      outfilename = arg
except GetoptError:
  print "bad command line arguments ... exiting ..."
  sys.exit(-1)
  
print "opening %s to read cflx" % infilename
nc = NC(infilename, 'r')
x = nc.variables["x"][:]
y = nc.variables["y"][:]
cflx = squeeze(nc.variables["cflx"][:])
print "  [cflx has max = %.2f and min = %.2f (m/a)]" % (cflx.max(),cflx.min())
nc.close()

mid = (len(x)-1)/2
plot(x[mid:]/1.e3,cflx[0,mid:],'k.-',markersize=10,linewidth=2)
hold(True)
plot(x[mid:]/1.e3,x[mid:] * 0.3,'r:',linewidth=1.5)
hold(False)
xlabel("y  ($\mathrm{km}$)",size=14)
ylabel(r"flux   ($\mathrm{m}^2\,\mathrm{a}^{-1}$)",size=14)
print "saving figure as %s ..." % outfilename
savefig(outfilename, dpi=300, facecolor='w', edgecolor='w')



