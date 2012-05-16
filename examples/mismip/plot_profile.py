#!/usr/bin/env python

import MISMIP

from pylab import figure, hold, plot, xlabel, ylabel, title, axis, vlines, savefig, text
from sys import exit

import numpy as np
import scipy
import scipy.optimize

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
                  help="Output image file name (e.g. -o foo.png)")
parser.add_option("-e", "--experiment", dest="experiment", type="string",
                  help="MISMIP experiment: 1a,1b,2a,2b,3a,3b (e.g. -e 1a)")
parser.add_option("-s", "--step", dest="step", type="int",
                  help="MISMIP step: 1,2,3,... (e.g. -s 1)")

(opts, args) = parser.parse_args()

infilename=""
outfilename=""

try:
    infilename = args[0]
except:
    print "ERROR: An input file is requied."
    exit(0)

tokens = infilename.split('_')
if tokens[0] == "ex":
    print "Using an -extra_file %s" % infilename
    tokens = tokens[1:]

try:
    model = tokens[0]
    experiment = tokens[1]
    mode = int(tokens[2][1])
    step = int(tokens[3][1])

    print "Plotting profile for model %s, experiment %s, grid mode %s, step %s" % (model, experiment, mode, step)
except:
    print "Can't parse file name %s" % infilename

    if opts.experiment is None:
        print "Please specify the experiment name (e.g. '-e 1a')."
        exit(0)
    else:
        experiment = opts.experiment

    if opts.step is None:
        print "Please specify the step (e.g. '-s 1')."
        exit(0)
    else:
        step = opts.step


if opts.output:
    outfilename=opts.output
else:
    import os.path
    outfilename=os.path.splitext(infilename)[0] + "-profile.pdf"

xg = MISMIP.x_g(experiment, step)

print "opening %s" % infilename
nc = NC(infilename)

x = nc.variables['x'][:]
thk = nc.variables['thk'][:]
usurf = nc.variables['usurf'][:]
topg = nc.variables['topg'][:]
mask = nc.variables['mask'][:]

if len(thk.shape) == 3:
    thk = thk[-1,1]
    usurf = usurf[-1,1]
    topg = topg[-1,1]
    mask = mask[-1,1]
else:
    thk = thk[1]
    usurf = usurf[1]
    topg = topg[1]
    mask = mask[1]

usurf = np.ma.array(usurf, mask=mask == 4)

lsrf = topg.copy()
lsrf[mask == 3] = -MISMIP.rho_i() / MISMIP.rho_w() * thk[mask == 3]
lsrf = np.ma.array(lsrf, mask=mask == 4)

def find_grounding_line(x, topg, thk, mask):
    # "positive" parts of x, topg, thk, mask
    x_p = x[x > 0]
    topg_p = topg[x > 0]
    thk_p  = thk[x > 0]
    mask_p = mask[x > 0]
    pass                                # unfinished

f = 1.0 - (- topg) / (thk * MISMIP.rho_i() / MISMIP.rho_w())

def func(x0):
    return scipy.interp(x0, x, f)

xg_PISM = scipy.optimize.fsolve(func, xg)

x /= 1e3
figure(1)
hold(True)
plot(x, np.zeros_like(x), ls='dotted', color='red')
plot(x, topg, color='black')
plot(x, usurf, 'o', color='blue', markersize=4)
plot(x, lsrf, 'o', color='blue', markersize=4)
xlabel('distance from the divide, km')
ylabel('elevation, m')
title("MISMIP experiment %s, step %d" % (experiment, step))
text(1200.0,4000.0,"$x_g$ (model) = %4.0f km" % (xg_PISM/1e3), color='r')
text(1200.0,3500.0,"$x_g$ (theory) = %4.0f km" % (xg/1e3))

_, _, ymin, ymax = axis(xmin=0, xmax=x.max())
vlines(xg/1e3, ymin, ymax, linestyles='dashed', color='black')
vlines(xg_PISM/1e3, ymin, ymax, linestyles='dashed', color='red')

print "saving figure as %s ..." % outfilename
savefig(outfilename)
