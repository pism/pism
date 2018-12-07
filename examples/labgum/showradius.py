#!/usr/bin/env python

# Copyright (C) 2013, 2014, 2016, 2017 the PISM Authors

# This script creates a graph of the modeled margin radius time series
# by using the ice_area_glacierized variable.  Compare Figure 4(a) in Sayag & Worster (2012).
# Try something like:
#    ./showradius.py -o foo.png ts_lab51.nc

import numpy as np
import pylab as plt
from sys import exit
try:
    import netCDF4 as netCDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

import argparse

parser = argparse.ArgumentParser(
    description='A script to show margin radius time series.  Requires one or more NetCDF files with "time" dimension and "time" and "ice_area_glacierized" variables.  Generates a .png image file with the graph.')
parser.add_argument('-o', '--outfile', metavar='FILENAME',
                    help='output file name (PNG)', default='foo.png')
parser.add_argument('-d', '--datafile', metavar='FILENAME',
                    help='data file name (ASCII with two columns:  t r_N)')
parser.add_argument('infiles', metavar='FILENAME', nargs='+',
                    help='input file name (NetCDF)')
args = parser.parse_args()

plt.figure(figsize=(12, 6))

for j in range(len(args.infiles)):
    nc = netCDF.Dataset(args.infiles[j], "r")
    t = nc.variables["time"][:]
    ice_area_glacierized = nc.variables["ice_area_glacierized"][:]
    nc.close()
    plt.loglog(t[t > 2], np.sqrt(ice_area_glacierized[t > 2] / np.pi) * 100.0, linewidth=2.0,  # after t=2s, and in cm
               label=args.infiles[j])

if args.datafile != None:
    A = np.loadtxt(args.datafile)
    data_t = A[:, 0]
    data_rN = A[:, 1]
    plt.loglog(data_t, 100.0 * data_rN, 'ko', label='observed', ms=4)  # cm versus s

plt.legend(loc='upper left')
plt.xticks([1.0, 10.0, 100.0, 1000.0])
plt.yticks([1.0, 10.0])
plt.axis([1.0, 1000.0, 1.0, 30.0])
plt.xlabel("t  (s)", size=14)
plt.ylabel(r"$r_N$  (cm)", size=14)
plt.grid(True)

plt.savefig(args.outfile)
