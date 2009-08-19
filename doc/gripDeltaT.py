#!/usr/bin/env python
"""This script creates gripDeltaT.png, the graph of the change of temperature from present, from the GRIP core.
The figure gripDeltaT.png is used in "Example: Modeling the Greenland ice sheet" in the User's Manual.
"""

from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from pylab import *

nc = NC("../examples/eisgreen/grip_dT.nc", "r")
t = nc.variables["t"][:]
delta_T = nc.variables["delta_T"][:]

figure(figsize=(12,6))
plot(t, delta_T, 'k')
xlabel("t (years before present)", size=16)
ylabel("$\Delta t$ (degrees C)", size=16)
grid(True)

savefig("figs/gripDeltaT.png")
