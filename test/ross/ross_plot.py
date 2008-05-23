#!/usr/bin/env python
# ross_plot.py plots the computed speed of Ross and shows observed values from RIGGS data. It is based on (and
# provided as an alternative to) plot_ross.m by ELB.
#
# Assumes that unnamed_diag.nc, containing the result of PISM Ross computation, and riggs_clean.dat, are is in the
# same directory.

# This script requires matplotlib for plotting (see http://matplotlib.sourceforge.net/ for details) and
# scipy.sandbox.delaunay for the 'natural neighbor' 2-d interpolation. This package is not built with scipy by
# default, but can be easily installed by downloading scipy sources, changing into the scipy/sandbox/delaunay
# directory and saying "python setup.py install". (This way it gets added to the global name-space; note the
# corresponding "import" line.)

# CK, May 22, 2008

from numpy import *
from pycdf import *
from pylab import *
from delaunay import *

# This is the name of the input file:
pism_output = 'unnamed_diag.nc'

# see 111by147.dat for these ranges
dlat = (-5.42445 - (-12.3325))/110;
gridlatext = linspace(-12.3325 - dlat * 46,-5.42445,147);
gridlon = linspace(-5.26168,3.72207,147);

[lat, lon] = meshgrid(gridlatext, gridlon)
glon = lon.flatten()
glat = lat.flatten()
tri = Triangulation(glon, glat)

# load RIGGS data FROM D. MACAYEAL TO ELB ON 19 DEC 2006.
RIGGS = load("riggs_clean.dat")

# read the PISM output
infile = CDF(pism_output)
H = infile.var("thk").get()[0]
mask = infile.var("mask").get()[0]
cbar = infile.var("cbar").get()[0]
ubar = infile.var("uvel").get()[0,:,:,0] * 31556926 # convert from meters per second to
vbar = infile.var("vvel").get()[0,:,:,0] * 31556926 # meters per year

# This is needed to only plot areas where H > 20 and mask == 0: 
cbar_masked = ma.array(cbar, mask = (mask == 1) + (H < 20))

# show computed speed as color
figure(1, figsize=(9,8));clf();hold(True)
pcolor(gridlon, gridlatext, cbar_masked); colorbar()

# compute grid lat and lon of RIGGS points (in deg,min,sec in .dat file); 
RIGGSlat = -(RIGGS[:,3] + RIGGS[:,4]/60 + RIGGS[:,5]/(60*60))
RIGGSlon = RIGGS[:,6] + RIGGS[:,7]/60 + RIGGS[:,8]/(60*60)
RIGGSlon = - RIGGSlon * RIGGS[:,9];  # RIGGS[:,9] is +1 if W, -1 if E

# throw out the ones which are not in model domain; 132 (135) remain
cRIGGS = tri.nn_interpolator(cbar.T.flatten())(RIGGSlon, RIGGSlat)
rig = RIGGS[cRIGGS > 0]; riglon = RIGGSlon[cRIGGS > 0]; riglat = RIGGSlat[cRIGGS > 0]

# add markers for RIGGS points, then quiver observed velocities
plot(riglon, riglat, '.k')
rigu = sin((pi/180)*rig[:,12]) * rig[:,10]
rigv = cos((pi/180)*rig[:,12]) * rig[:,10]
quiver(riglon, riglat, rigu, rigv, color='black')

# quiver the computed velocities at the same points; note reversal of u,v in model
uATrig = tri.nn_interpolator(vbar.T.flatten())(riglon, riglat)
vATrig = tri.nn_interpolator(ubar.T.flatten())(riglon, riglat)
quiver(riglon, riglat, uATrig, vATrig, color='red')
axis([-5.26168, 3.72207, -13, -5.42445])
xlabel('RIGGS grid longitude (deg E)'); ylabel('RIGGS grid latitude (deg N)')
title('Color is speed in m/a.\n Arrows are observed (black) and computed (red) velocities at RIGGS points.')
show()
#savefig("rossquiver.pdf")

figure(2);clf();hold(True)
# report results comparable to Table 1 in (MacAyeal et al 1996)
ChiSqrActual = sum( ((uATrig - rigu)**2 + (vATrig - rigv)**2) / (30**2) )
print "ChiSqr = %f" % (ChiSqrActual * (156.0/132.0))
print "Maximal computer speed is %f." % (cbar.max())

# show observed versus computed scatter plot as in Figure 2 in (MacAyeal et al 1996)
plot(sqrt(uATrig**2 + vATrig**2), sqrt(rigu**2 + rigv**2), '.k')
plot([0, 1000],[0, 1000], 'k')
xlabel('PISM computed speed (m/a)'); ylabel('RIGGS observed speed (m/a)')
show()
#savefig("rossscatter.pdf")
