#!/usr/bin/env python

from pylab import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
import sys

from optparse import OptionParser

parser = OptionParser()
parser.usage = "%prog <input file> [options]"
parser.description = "Plots the ice flux as a function of the distance from the divide."
parser.add_option("-o", "--output", dest="output", type="string",
                  help="Output file name")

(opts, args) = parser.parse_args()

infilename=""
outfilename=""

try:
    infilename = args[0]
except:
    print "ERROR: An input file is requied."
    exit(0)

if opts.output:
    outfilename=opts.output
else:
    import os.path
    outfilename=os.path.splitext(infilename)[0] + "-flux.pdf"

print "opening %s to read cflx" % infilename
nc = NC(infilename, 'r')
x = nc.variables["x"][:]
cflx = squeeze(nc.variables["cflx"][:])
nc.close()

# extract the last record from an -extra_file, if necessary:
if len(cflx.shape) == 3:
    cflx = cflx[-1]

print "  [cflx has max = %.2f and min = %.2f (m/a)]" % (cflx.max(),cflx.min())

mid = (len(x)-1)/2
plot(x[mid:]/1.e3,cflx[0,mid:],'k.-',markersize=10,linewidth=2)  # FIXME: variable order matters here

hold(True)
plot(x[mid:]/1.e3,x[mid:] * 0.3,'r:',linewidth=1.5)

xlabel("x  ($\mathrm{km}$)",size=14)
ylabel(r"flux   ($\mathrm{m}^2\,\mathrm{a}^{-1}$)",size=14)

print "saving figure as %s ..." % outfilename
savefig(outfilename, dpi=300, facecolor='w', edgecolor='w')

