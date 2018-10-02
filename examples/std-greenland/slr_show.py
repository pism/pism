#!/usr/bin/env python

# Copyright (C) 2011-2012, 2014, 2016, 2017 The PISM Authors

# script to generate figure: results from SeaRISE experiments
# usage:  if UAFX_G_D3_C?_??.nc are result NetCDF files then do
# $ slr_show.py -m UAFX

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

from numpy import zeros
import pylab as plt
from optparse import OptionParser

parser = OptionParser()
parser.usage = "usage: %prog [options]"
parser.description = "A script for PISM output files to show time series plots using pylab."
parser.add_option("-a", dest="t_a", type="int",
                  help="start year, in years since 2004, default = 0", default=0)
parser.add_option("-e", dest="t_e", type="int",
                  help="end year, in years since 2004, default = 500", default=500)
parser.add_option("-m", "--model", dest="model",
                  help="choose experiment, default UAF1", default="UAF1")


(options, args) = parser.parse_args()
model = options.model
t_a = options.t_a
t_e = options.t_e

# first name in this list is CONTROL
NCNAMES = [model + "_G_D3_C1_E0.nc", model + "_G_D3_C2_E0.nc", model + "_G_D3_C3_E0.nc", model + "_G_D3_C4_E0.nc", model + "_G_D3_C1_S1.nc", model +
           "_G_D3_C1_S2.nc", model + "_G_D3_C1_S3.nc", model + "_G_D3_C1_M1.nc", model + "_G_D3_C1_M2.nc", model + "_G_D3_C1_M3.nc", model + "_G_D3_C1_T1.nc"]

# labels
labels = ["AR4 A1B", "AR4 A1B 1.5x", "AR4 A1B 2x", "2x basal sliding", "2.5x basal sliding",
          "3x basal sliding", "2 m/a bmr", "20 m/a bmr", "200 m/a bmr", "AR4 A1B + 2x sliding"]
# line colors
colors = ['#984EA3',  # violet
          '#984EA3',  # violet
          '#984EA3',  # violet
          '#FF7F00',  # orange
          '#FF7F00',  # orange
          '#FF7F00',  # orange
          '#377EB8',  # light blue
          '#377EB8',  # light blue
          '#377EB8',  # light blue
          '#4DAF4A']  # green

dashes = ['-', '--', '-.', '-', '--', '-.', '-', '--', '-.', '-']

print("control run name is " + NCNAMES[0])

n = len(NCNAMES)
nc0 = CDF(NCNAMES[0], 'r')
try:
    t_units = nc0.variables['tseries'].units
    t = nc0.variables['tseries'][t_a:t_e]
except:
    t_units = nc0.variables['time'].units
    t = nc0.variables['time'][t_a:t_e]
nc0.close()

# convert to years if t is in seconds
if (t_units.split()[0] == ('seconds' or 's')):
    t /= 3.15569259747e7

ice_volume_glacierized = zeros((len(t), n))
ivolshift = zeros((len(t), n - 1))

for j in range(n):
    nc = CDF(NCNAMES[j], 'r')
    ice_volume_glacierized[:, j] = nc.variables['ice_volume_glacierized'][t_a:t_e]
    nc.close()

for j in range(n - 1):
    ivolshift[:, j] = ice_volume_glacierized[:, j + 1] - ice_volume_glacierized[:, 0]

# "2,850,000 km3 of ice were to melt, global sea levels would rise 7.2 m"
scale = 7.2 / 2.850e6


# screen plot with high contrast
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='0.15')
for j in range(n - 1):
    ax.plot(t, -(ivolshift[:, j] / 1.0e9) * scale, dashes[j], color=colors[j], linewidth=3)
ax.set_xlabel('years from 2004')
ax.set_ylabel('sea level rise relative to control (m)')
ax.legend(labels, loc='upper left')
ax.grid(True, color='w')

plt.show()


# line colors
colors = ['#984EA3',  # violet
          '#984EA3',  # violet
          '#984EA3',  # violet
          '#FF7F00',  # orange
          '#FF7F00',  # orange
          '#FF7F00',  # orange
          '#084594',  # dark blue
          '#084594',  # dark blue
          '#084594',  # dark blue
          '#4DAF4A']  # green

# print plot with white background
fig = plt.figure()
ax = fig.add_subplot(111)
for j in range(n - 1):
    ax.plot(t, -(ivolshift[:, j] / 1.0e9) * scale, dashes[j], color=colors[j], linewidth=2)
ax.set_xlabel('years from 2004')
ax.set_ylabel('sea level rise relative to control (m)')
ax.legend(labels, loc='upper left')
ax.grid(True)
plt.savefig(model + '_slr.pdf')
