#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from argparse import ArgumentParser

# Import all necessary modules here so that if it fails, it fails early.
try:
    import netCDF4 as NC
except:
    import netCDF3 as NC

# Set up the option parser
parser = ArgumentParser()
parser.description = "A script to plot results of the ROSS example."
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE

if len(args) == 1:
    pism_output = args[0]
else:
    print("wrong number of arguments, 1 expected, %i given" % int(len(args)))
    import sys
    exit(1)
    
try:
    nc = NC.Dataset(pism_output, 'r')
except:
    print("file %s not found" % pism_output)
    import sys
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
plt.savefig('rossquiver.png', dpi=300)

# do the scatter plot
magnitude = np.sqrt(np.abs(u[::s,::s])**2 + np.abs(v[::s,::s])**2)
bc_magnitude = np.sqrt(np.abs(u_bc[::s,::s])**2 + np.abs(v_bc[::s,::s])**2)

max_velocity = np.maximum(magnitude.max(), bc_magnitude.max())
plt.figure(2)
plt.hold(True)
plt.scatter(magnitude, bc_magnitude, color='black')
plt.plot([0, max_velocity], [0, max_velocity], color='black', ls='--')
plt.axis(xmin=0, xmax=max_velocity, ymin=0, ymax=max_velocity)
plt.xlabel('modeled speed')
plt.ylabel('observed speed')
plt.title("Observed versus modeled speed (m/year), at points in quiver plot")
print "saving figure 'rossscatter.png'"
plt.savefig('rossscatter.png', dpi=300)
