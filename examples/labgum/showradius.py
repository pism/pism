#!/usr/bin/env python

# Copyright (C) 2013 the PISM Authors

# This script creates a graph of the modeled margin radius time series
# by using the iarea variable.  Compare Figure 4(a) in Sayag & Worster (2012).
# Try something like:
#    ./showradius.py -o foo.png ts_lab51.nc

import numpy as np
import pylab as plt
from sys import exit
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF

import argparse

parser = argparse.ArgumentParser(description='A script to show margin radius time series.  Requires a NetCDF file with "time" dimension and "time" and "iarea" variables.  Generates a .png image file with the graph.')
parser.add_argument('-o', '--outfile', metavar='FILENAME',
                   help='output file name (PNG)',default='foo.png')
parser.add_argument('-d', '--datafile', metavar='FILENAME',
                   help='data file name (ASCII with two columns:  t r_N)')
parser.add_argument('infile', metavar='FILENAME',
                   help='input file name (NetCDF)')
args = parser.parse_args()

nc = netCDF.Dataset(args.infile, "r")
t = nc.variables["time"][:]
iarea = nc.variables["iarea"][:]
nc.close()

plt.figure(figsize=(12,6))
plt.loglog(t[t>2], np.sqrt(iarea[t>2]/np.pi) * 100.0, 'b', linewidth=2.0)  # after t=2s, and in cm
#plt.legend(('Q_0','Regression'), loc='upper left')

if args.datafile != None:
  A = np.loadtxt(args.datafile)
  data_t = A[:,0]
  data_rN = A[:,1]
  plt.hold(True)
  plt.loglog(data_t, 100.0 * data_rN, 'rx')  # cm versus s
  plt.hold(False)

plt.xticks([1.0, 10.0, 100.0, 1000.0])
plt.yticks([1.0, 10.0])
plt.axis([1.0, 1000.0, 1.0, 30.0])
plt.xlabel("t  (s)", size=14)
plt.ylabel(r"$r_N$  (cm)", size=14)
plt.grid(True)

plt.savefig(args.outfile)

