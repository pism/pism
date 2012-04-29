#!/usr/bin/env python
# ross_plot.py plots the computed Ross ice shelf speed compared to observed
# values from RIGGS data.
#
# This script depends on the following Python packages:
#   1) numpy (see http://numpy.scipy.org/; python-numpy is Debian package)
#   2) matplotlib (see http://matplotlib.sourceforge.net/; python-matplotlib 
#        is Debian package)
#   3) netcdf4-python (see http://code.google.com/p/netcdf4-python/; no Debian package)
#
# CK 27may08, ..., 12jan10
# ELB 15feb11

from numpy import ma, loadtxt, squeeze, shape, reshape, linspace, tile, repeat, sin, pi, cos, sqrt, maximum
from pylab import figure, clf, hold, pcolor, colorbar, plot, quiver, axis, xlabel, ylabel, savefig, show, title
from getopt import getopt, GetoptError
from sys import argv, exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    from matplotlib.delaunay import Triangulation
except:
    print "ERROR: matplotlib.delaunay  not installed (?)"
    print "Please upgrade matplotlib"
    exit(1)

seconds_per_year = 3.1556926e7

# process command line arguments
try:
    opts, args = getopt(argv[1:], "", ["pism-output=", "riggs="])
    # defaults:
    pism_output = "rossComputed.nc"
    riggs_file = "riggs_clean.dat"
    for opt, arg in opts:
        if opt in ("--pism-output"):
            pism_output = arg
        if opt in ("--riggs"):
            riggs_file = arg
except GetoptError:
    print """
Options:
   --pism-output=<PISM output .nc file>:  specifies the NetCDF file with PISM output
   --riggs=<RIGGS data file>: specifies the data file with RIGGS points
"""
    exit(-1)

# load RIGGS data FROM D. MACAYEAL TO ELB ON 19 DEC 2006.
try:
    print "Loading RIGGS data from '%s'..." % (riggs_file),
    RIGGS = loadtxt(riggs_file)  # pylab now suggests numpy.loadtxt instead of 
                                 #    pylab's "load"
    print "done."
except IOError:
    print """ERROR!\nMake sure that '%s' is in the expected location
     and try again.  Exiting...""" % (riggs_file)
    exit(-1)

# load the PISM output
try:
    print "Loading PISM output from '%s'..." % (pism_output),
    infile = NC(pism_output, 'r')
except Exception:
    print """ERROR!\nSpecify NetCDF file from PISM run with -p.
    See ross_plot.py --help and User's Manual.  Exiting..."""
    exit(-1)

def get_var(nc, name):
    """Get a 2D field from a NetCDF file, transposing if necessary. The
    returned array will always have the (y,x) variable order."""

    var = nc.variables[name]
    dims = var.dimensions
    if dims.index('x') < dims.index('y'):
        # transpose
        return squeeze(var[:]).T
    else:
        return squeeze(var[:])

H         = get_var(infile, "thk")
mask      = get_var(infile, "mask")
cbar      = get_var(infile, "cbar")
ubar      = get_var(infile, "u_ssa")
vbar      = get_var(infile, "v_ssa")
pismaccur = get_var(infile, "accur")
print "done."

# this may not work if axis order is flipped; works r1443 with Pism_FAST_WRITE=2
(Mx,My) = shape(cbar)

# see 111by147.dat for these ranges
dlat = (-5.42445 - (-12.3325))/110
gridlatext = linspace(-12.3325 - dlat * 46,-5.42445,Mx)
gridlon = linspace(-5.26168,3.72207,My)

# need RIGGS lat,lon in different forms
glat = repeat(gridlatext, My)
glon = tile(gridlon, Mx)
reglat = reshape(glat,(Mx,My))
reglon = reshape(glon,(Mx,My))

# Plot areas where thickness and PISM mask are appropriate,
# and filter out RIGGS points that are outside the model domain,
# and use location of calving front from pismaccur flag, which is
# interpolated by PISM from original ross.nc file "accur" flag.
cbar_masked = ma.array(cbar, mask = (mask != 3) + (H < 20.0)
             + ((pismaccur < 0.01) & (reglat < -11.0) & (reglon < 1.0) ) )

# show computed speed as color
figure(1, figsize=(9,8));clf();hold(True)
pcolor(gridlon, gridlatext, cbar_masked); colorbar()

# compute grid lat and lon of RIGGS points (in deg,min,sec in .dat file); 
RIGGSlat = -(RIGGS[:,3] + RIGGS[:,4]/60 + RIGGS[:,5]/(60*60))
RIGGSlon = RIGGS[:,6] + RIGGS[:,7]/60 + RIGGS[:,8]/(60*60)
RIGGSlon = - RIGGSlon * RIGGS[:,9];  # RIGGS[:,9] is +1 if W, -1 if E

# throw out the ones which are not in model domain; 132 (131?) remain
cbar_masked = cbar_masked.filled(-20)

# triangulate data
tri = Triangulation(glon, glat)
cRIGGS = tri.nn_interpolator(cbar_masked.flat)(RIGGSlon, RIGGSlat)
rig = RIGGS[cRIGGS > 0]
riglon = RIGGSlon[cRIGGS > 0]
riglat = RIGGSlat[cRIGGS > 0]

# add markers for RIGGS points, then quiver observed velocities
plot(riglon, riglat, '.k')
rigu = sin((pi/180)*rig[:,12]) * rig[:,10]
rigv = cos((pi/180)*rig[:,12]) * rig[:,10]
quiver(riglon, riglat, rigu, rigv, color='black')

# quiver the computed velocities at the same points
uATrig = tri.nn_interpolator(ubar.flat)(riglon, riglat)
vATrig = tri.nn_interpolator(vbar.flat)(riglon, riglat)
quiver(riglon, riglat, uATrig, vATrig, color='red')
axis([-5.26168, 3.72207, -12.75, -5.42445])
xlabel('RIGGS grid longitude (deg E)'); ylabel('RIGGS grid latitude (deg N)')
title("""Color is speed in m/a.\n Arrows are observed (black) and computed 
(red) velocities at RIGGS points.""")

# to report results comparable to Table 1 in (MacAyeal et al 1996)
#ChiSqrActual = sum( ((uATrig - rigu)**2 + (vATrig - rigv)**2) / (30**2) )
#print "chi^2 = %f" % (ChiSqrActual * (156.0/132.0))
#print "maximum computed ice shelf speed is  %f." % (cbar.max())

# show observed versus computed scatter plot as in Figure 2 in (MacAyeal et al 1996)
figure(2);clf();hold(True)
pism_result = sqrt(uATrig**2 + vATrig**2)
riggs_data = sqrt(rigu**2 + rigv**2)
max_speed = maximum(pism_result.max(), riggs_data.max())

plot(pism_result, riggs_data, '.k')
plot([0, max_speed],[0, max_speed], 'k')
axis(xmax=max_speed, ymax=max_speed)
xlabel('PISM computed speed (m/a)'); ylabel('RIGGS observed speed (m/a)')

print "pausing to show figures ..."
show()

