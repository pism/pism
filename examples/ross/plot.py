#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from argparse import ArgumentParser

# Import all necessary modules here so that if it fails, it fails early.
try:
    import netCDF4 as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

# Set up the option parser
parser = ArgumentParser()
parser.description = "A script to plot results of the ROSS example."
parser.add_argument("FILE", nargs='*')

options = parser.parse_args()
args = options.FILE

if len(args) == 1:
    pism_output = args[0]
    plotname = pism_output.split('.')[0]
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
    return np.ma.array(a, mask=mask != 3)


# load data and restrict velocities to floating areas
seconds_per_year = 3.1556926e7
x = get('x')
y = get('y')
velsurf_mag = get('velsurf_mag')
mask = get('mask')
u = floating(get('u_ssa'))
v = floating(get('v_ssa'))
# B.C.s are observations, so a PISM output file contains everything we need
u_bc = floating(get('u_ssa_bc'))
v_bc = floating(get('v_ssa_bc'))

plt.clf()

f, (a0, a1) = plt.subplots(1, 2,
                           gridspec_kw={'width_ratios': [1.2, 1]},
                           figsize=(16, 8))

# mark the grounding line
a0.contour(x, y, mask, [2.5], colors="black", lw=2)

# plot velsurf_mag using log color scale
p = a0.pcolormesh(x, y, velsurf_mag, norm=colors.LogNorm(vmin=1, vmax=1.5e3))
# add a colorbar:
f.colorbar(p, ax=a0, extend='both', ticks=[1, 10, 100, 500, 1000], format="%d")

# quiver plot of velocities
s = 10                                  # stride in grid points
a0.quiver(x[::s], y[::s], u_bc[::s, ::s], v_bc[::s, ::s], color='white')
a0.quiver(x[::s], y[::s], u[::s, ::s], v[::s, ::s], color='black')
a0.set_xticks([])
a0.set_yticks([])
a0.set_title("Ross ice velocity (m/year)\nwhite=observed, black=model")

# do the scatter plot
magnitude = np.sqrt(np.abs(u[::s, ::s]) ** 2 + np.abs(v[::s, ::s]) ** 2)
bc_magnitude = np.sqrt(np.abs(u_bc[::s, ::s]) ** 2 + np.abs(v_bc[::s, ::s]) ** 2)

max_velocity = np.maximum(magnitude.max(), bc_magnitude.max())

a1.scatter(magnitude, bc_magnitude, marker=".")
a1.plot([0, max_velocity], [0, max_velocity], color='black', ls='--')
a1.axis(xmin=0, xmax=max_velocity, ymin=0, ymax=max_velocity)
a1.set_xlabel('modeled speed')
a1.set_ylabel('observed speed')
a1.set_title("Observed versus modeled speed (m/year)\nat points in quiver plot")

output_filename = plotname + '.png'

f.tight_layout()

f.savefig(output_filename, dpi=300)

print("saving figure %s" % output_filename)
