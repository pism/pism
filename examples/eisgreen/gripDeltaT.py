#!/usr/bin/env python
"""This script creates gripDeltaT.pdf, the graph of the change of temperature
from present, from the GRIP core.  The figure gripDeltaT.pdf is used in
"Example: Modeling the Greenland ice sheet" in the User's Manual.
"""

from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from pylab import *

nc = NC("grip_dT.nc", "r")
t = nc.variables["t"][:]
delta_T = nc.variables["delta_T"][:]

fig_width  = 6.68  # width in inches
fig_height = 2.00  # height in inches
fig_size   = [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size}
rcParams.update(params)


figure()
plot(-t, delta_T, 'k')
xlabel("t (years before present)")
ylabel("$\Delta T$ (degrees C)")
grid(True)

savefig("gripDeltaT.pdf",bbox_inches='tight')
