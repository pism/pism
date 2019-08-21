#!/usr/bin/env python

from PISMNC import PISMDataset as PNC
import numpy as np
from sys import exit

# Import all necessary modules here so that if it fails, it fails early.
try:
    import netCDF4 as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

inname = "pismnbreen.nc"
outname = "fakesummerevent.nc"

try:
    innc = NC.Dataset(inname, 'r')
except:
    print("file %s not found" % inname)
    exit(1)

try:
    nc = PNC(outname, 'w', format='NETCDF3_CLASSIC')
except:
    print("can't open file %s for writing" % outname)
    exit(1)


def get(name):
    global innc
    return np.squeeze(innc.variables[name][:])


x = get('x')
y = get('y')
bmelt = get('basal_melt_rate_grounded')
Mx = len(x)
My = len(y)
zero = np.zeros((My, Mx))

nc.create_dimensions(x, y, time_dependent=True, use_time_bounds=True)


def drainage(t):
    """time-dependence of bogus summer runoff event in m/a: a positive wavepacket"""
    return np.exp(-(t - 180.0) ** 2 / 80.0) * 20.0 * (np.cos(0.2 * t * 2 * 3.14159) + 1.0)


year = 2012
nc.variables['time'].units = "days since %d-1-1" % (year)

# generate space-time bogus summer runoff event; mask where bmelt > 0
for a in range(1, 366):
    nc.append_time(a, (a, a + 1))
    inputthisday = (zero + drainage(np.double(a))) * (bmelt > 0)
    nc.write("inputtobed", inputthisday, True)

# Set attributes
inputtobed = nc.variables["inputtobed"]
inputtobed.units = "m / year"

nc.close()
innc.close()
