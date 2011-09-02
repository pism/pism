#!/usr/bin/env python

from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from pylab import *

nc_no = NC("ts_no_force.nc", "r")
nc_default = NC("ts_with_force.nc", "r")
nc_weak = NC("ts_weak_force.nc", "r")
nc_strong = NC("ts_strong_force.nc", "r")

stride = 10

try:
    t = nc_no.variables["time"][::stride]
except:
    t = nc_no.variables["t"][::stride]
    
ivol_no = nc_no.variables["ivol"][::stride]
ivol_default = nc_default.variables["ivol"][::stride]
ivol_weak = nc_weak.variables["ivol"][::stride]
ivol_strong = nc_strong.variables["ivol"][::stride]

figure()
plot(t, ivol_no * 1.0e-15, '-o',
     t, ivol_default * 1.0e-15, '-d',
     t, ivol_weak * 1.0e-15, '-s',
     t, ivol_strong * 1.0e-15, '-^',
     t, 2.82528 * ones(shape(t)),'--k')
legend( ('no forcing',
         'alpha = 0.01 (default)',
         'alpha = 0.005',
         'alpha = 0.05',
         'target volume'),
         loc='right')
xlabel("t (years)", size=16)
ylabel("ice volume (10$^6$ km$^3$)", size=16)
grid(True)

savefig("ivol_force_to_thk.png", dpi=100)

figure()
diffusivity_no = nc_no.variables["max_diffusivity"][::stride]
diffusivity_default = nc_default.variables["max_diffusivity"][::stride]
diffusivity_weak = nc_weak.variables["max_diffusivity"][::stride]
diffusivity_strong = nc_weak.variables["max_diffusivity"][::stride]

plot(t, diffusivity_no, '-o',
     t, diffusivity_default, '-d',
     t, diffusivity_weak, '-s',
     t, diffusivity_strong, '-^')
legend( ('no forcing',
         'alpha = 0.01 (default)',
         'alpha = 0.005',
         'alpha = 0.05',
         'target volume'),
         loc='right')
xlabel("t (years)", size=16)
ylabel("maximum diffusivity (m$^{2}$ s$^{-1}$)", size=16)
grid(True)

savefig("diffusivity_force_to_thk.png", dpi=100)

