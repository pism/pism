#!/usr/bin/env python

from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from pylab import *

nc_no = NC("ts_no_force.nc", "r")
nc_with = NC("ts_with_force.nc", "r")
nc_weak = NC("ts_weak_force.nc", "r")

t = nc_no.variables["t"][:]

ivol_no = nc_no.variables["ivol"][:]
ivol_with = nc_with.variables["ivol"][:]
ivol_weak = nc_weak.variables["ivol"][:]

plot(t, ivol_no * 1.0e-15, '-o',
     t, ivol_with * 1.0e-15, '-o',
     t, ivol_weak * 1.0e-15, '-o',
     t, 2.82528 * ones(shape(t)),'--k')
legend( ('no forcing',
         'alpha = 0.002 (default)',
         'alpha = 0.0002',
         'target volume'),
         loc='right')
xlabel("t (years)", size=16)
ylabel("ice volume (10^6 km^3)", size=16)
grid(True)

savefig("ivol_force_to_thk.png", dpi=100)

