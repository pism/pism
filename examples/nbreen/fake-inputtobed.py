#!/usr/bin/env python3

from PISMNC import PISMDataset as PNC
import numpy as np
import netCDF4 as NC

inname = "pismnbreen.nc"
outname = "fakesummerevent.nc"

innc = NC.Dataset(inname, 'r')

nc = PNC(outname, 'w', format='NETCDF3_CLASSIC')

def get(name):
    global innc
    return np.squeeze(innc.variables[name][:])

x = get('x')
y = get('y')
bmelt = get('basal_melt_rate_grounded')
zero = np.zeros_like(bmelt)

nc.create_dimensions(x, y, time_dependent=True, use_time_bounds=True)

def drainage(t):
    """time-dependence of bogus summer runoff event in m/a: a positive wavepacket"""
    C = np.cos(0.2 * t * 2 * 3.14159) + 1.0
    return np.exp(-(t - 180.0) ** 2 / 80.0) * 20.0 * C

year = 2012
nc.variables['time'].units = "days since {}-1-1".format(year)

water_density = 1000.0
# generate space-time bogus summer runoff event; mask where bmelt > 0
for a in range(1, 366):
    nc.append_time(a, (a, a + 1))
    inputthisday = (zero + drainage(np.double(a))) * (bmelt > 0) * water_density
    nc.write("water_input_rate", inputthisday, True)

# Set attributes
nc.variables["water_input_rate"].units = "kg m^-2 year^-1"

nc.close()
innc.close()
