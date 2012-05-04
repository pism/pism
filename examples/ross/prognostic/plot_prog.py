#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from getopt import getopt, GetoptError
from sys import argv, exit

# Import all necessary modules here so that if it fails, it fails early.
try:
    import netCDF4 as NC
except:
    import netCDF3 as NC

# process command line arguments
try:
    opts, args = getopt(argv[1:], "", ["pism-output="])
    # defaults:
    pism_output = "out.nc"
    for opt, arg in opts:
        if opt in ("--pism-output"):
            pism_output = arg
except GetoptError:
    print """
Options:
   --pism-output=<PISM output .nc file>:  specifies the NetCDF file with PISM output
"""
    exit(-1)


try:
    nc = NC.Dataset(pism_output, 'r')
except:
    print "file %s missing; use --pism-output=foo.nc to read from foo.nc" % pism_output
    exit(1)


def get(name):
    global nc
    return np.squeeze(nc.variables[name][:])

def floating(a):
    global mask
    return np.ma.array(a, mask=mask!=3)

# load data and restrict velocities to floating areas
seconds_per_year = 3.1556926e7
x     = get('x')
y     = get('y')
csurf = get('csurf')
mask  = get('mask')
u     = floating(get('u_ssa'))
v     = floating(get('v_ssa'))
# B.C.s are observations, so a PISM output file contains everything we need
u_bc  = floating(get('u_ssa_bc')) * seconds_per_year # convert to m/year
v_bc  = floating(get('v_ssa_bc')) * seconds_per_year

# dark background
plt.figure(1)
plt.subplot(111, axisbg='0.5')

# mark the grounding line
plt.contour(x, y, mask, [2.5], colors="black", lw=2)

# plot csurf using log color scale
plt.pcolormesh(x, y, csurf, norm=colors.LogNorm(vmin=1, vmax=1.5e3))
# add a colorbar:
plt.colorbar(extend='both', ticks=[1, 10, 100, 250, 500, 1000], format="%d")

# quiver plot of velocities
s = 10                                  # stride in grid points
plt.quiver(x[::s], y[::s], u_bc[::s,::s], v_bc[::s,::s], color='white')
plt.quiver(x[::s], y[::s], u[::s,::s], v[::s,::s], color='black')
plt.xticks([])
plt.yticks([])
plt.title(r"Ross ice velocity (m/year); white=observed, black=model")
print "saving figure 'rossquiver.png'"
plt.savefig('prog_rossquiver.png', dpi=300)

# do the scatter plot
magnitude = np.sqrt(np.abs(u[::s,::s])**2 + np.abs(v[::s,::s])**2)
bc_magnitude = np.sqrt(np.abs(u_bc[::s,::s])**2 + np.abs(v_bc[::s,::s])**2)

magnitude_no0 =np.ma.array(magnitude, mask = (magnitude==0.0))
bc_magnitude_no0 =np.ma.array(bc_magnitude, mask = (bc_magnitude==0.0))

max_velocity = np.maximum(magnitude_no0.max(), bc_magnitude_no0.max())
plt.figure(2)
plt.hold(True)
plt.scatter(magnitude_no0, bc_magnitude_no0, color='black')
plt.plot([0, max_velocity], [0, max_velocity], color='black', ls='--')
plt.axis(xmin=0, xmax=max_velocity, ymin=0, ymax=max_velocity)
plt.xlabel('modeled speed')
plt.ylabel('observed speed')
plt.title("Observed versus modeled speed (m/year), at points in quiver plot")
print "saving figure 'rossscatter.png'"
plt.savefig('prog_rossscatter.png', dpi=300)
