#!/usr/bin/env python
"""This script creates ivolboth.png, the graph of the modeled volume time series
for both 20km and 10km grid calculations for SeaRISE-Greenland.  This figure
appears in the "Getting Started" section of the User's Manual.
"""

from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from pylab import *

# generate "whole" time series this way:
#   $ ncrcat ts_g20km_m40ka.nc ts_g20km_m30ka.nc ts_g20km_m20ka.nc \
#       ts_g20km_m10ka.nc ts_g20km_0.nc -o ts_g20km_whole.nc

nc20km = NC("ts_g20km_whole.nc", "r")
t20km = nc20km.variables["t"][:]
ivol20km = nc20km.variables["ivol"][:]
nc20km.close()

nc10km = NC("ts_g10km_whole.nc", "r")
t10km = nc10km.variables["t"][:]
ivol10km = nc10km.variables["ivol"][:]
nc10km.close()

vfactor = 1.0e6 * 1.0e9

figure(figsize=(12,6))
plot(t20km, ivol20km / vfactor, 'b', linewidth=1.0), hold(True)
plot(t10km, ivol10km / vfactor, 'r', linewidth=2.0), hold(False)
legend(('20 km','10 km'))
xlabel("t (years)", size=16)
ylabel("volume (10^6 km^3)", size=16)
grid(True)

last=-10000.0
t20km_last = t20km[t20km >= last]
ivol20km_last = ivol20km[t20km >= last]
t10km_last = t10km[t10km >= last]
ivol10km_last = ivol10km[t10km >= last]

axesinset = axes([0.55, 0.20, 0.33, 0.25],axisbg='w')
plot(t20km_last,ivol20km_last / vfactor, 'b', linewidth=1.0), hold(True)
plot(t10km_last,ivol10km_last / vfactor, 'r', linewidth=2.0), hold(False)
setp(axesinset)
#setp(axesinset,xticks=arange(95.e3,101.e3,1.e3),
#     xticklabels=('95','96','97','98','99','100'))

show()

#savefig("ivolboth.png")

