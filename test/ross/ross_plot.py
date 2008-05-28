#!/usr/bin/env python
# ross_plot.py plots the computed Ross ice shelf speed compared to  observed 
# values from RIGGS data.
#
# This script depends on the following Python packages:
#   1) numpy (see http://numpy.scipy.org/; python-numpy is Debian package)
#   2) matplotlib (see http://matplotlib.sourceforge.net/; python-matplotlib 
#        is Debian package)
#   3) pycdf (see http://pysclint.sourceforge.net/pycdf/; no Debian package)
#   4) scikits.delaunay (for data regridding; see 
#        http://scipy.org/scipy/scikits for more information and to learn 
#        what Scipy Toolkits are.)
#
# To install scikits.delaunay, first install setuptools (python-setuptools 
# is Debian package).  Then do
#
#    $ svn co http://svn.scipy.org/svn/scikits/trunk/delaunay
#    $ cd delaunay; sudo python setup.py install
#   
# Note there is an alternative Matlab script ross_plot.m
#
# CK 27may08

from numpy import *
from pylab import *
from getopt import getopt, GetoptError
from sys import argv, exit
from pycdf import CDF, CDFError
from scikits.delaunay import *

# process command line arguments
try:
    opts, args = getopt(argv[1:], "p:r:", ["pism-output=", "riggs="])
    # defaults:
    pism_output = "unnamed_diag.nc"
    riggs_file = "riggs_clean.dat"
    for opt, arg in opts:
        if opt in ("-p", "--pism-output"):
            pism_output = arg
        if opt in ("-r", "--riggs"):
            riggs_file = arg
except GetoptError:
    print """Options:
               --pism-output=<PISM output file> or -p <PISM output file>
                 specifies the NetCDF file with PISM output

               --riggs=<RIGGS data file> or -r <RIGGS data file>
                 specifies the data file with RIGGS points."""
    exit(-1)

# load RIGGS data FROM D. MACAYEAL TO ELB ON 19 DEC 2006.
try:
    print "Loading RIGGS data from '%s'..." % (riggs_file),
    RIGGS = load(riggs_file)
    print "done."
except IOError:
    print """ERROR!\nMake sure that '%s' is in the expected location
     and try again.  Exiting...""" % (riggs_file)
    exit(-1)

# load the PISM output
try:
    print "Loading PISM output from '%s'..." % (pism_output),
    infile = CDF(pism_output)
    H = infile.var("thk").get()[0]
    mask = infile.var("mask").get()[0]
    cbar = infile.var("cbar").get()[0]
    ubar = infile.var("uvel").get()[0,:,:,0] * (60*60*24*365) # convert from 
    vbar = infile.var("vvel").get()[0,:,:,0] * (60*60*24*365) # m s^-1 to m a^-1
    print "done."
except CDFError:
    print """ERROR!\nSpecify NetCDF file from PISM run with -p.
    See ross_plot.py --help and User's Manual.  Exiting..."""
    exit(-1)

# see 111by147.dat for these ranges
dlat = (-5.42445 - (-12.3325))/110;
gridlatext = linspace(-12.3325 - dlat * 46,-5.42445,147);
gridlon = linspace(-5.26168,3.72207,147);

# triangulate data
glon = tile(gridlon, 147); glat = repeat(gridlatext, 147)
tri = Triangulation(glon, glat)

# This is needed to only plot areas where H >= 20 and mask == 0
# and filter out RIGGS points that are outside the model domain.
cbar_masked = ma.array(cbar, mask = (mask == 1) + (H < 20))

# show computed speed as color
figure(1, figsize=(9,8));clf();hold(True)
pcolor(gridlon, gridlatext, cbar_masked); colorbar()

# compute grid lat and lon of RIGGS points (in deg,min,sec in .dat file); 
RIGGSlat = -(RIGGS[:,3] + RIGGS[:,4]/60 + RIGGS[:,5]/(60*60))
RIGGSlon = RIGGS[:,6] + RIGGS[:,7]/60 + RIGGS[:,8]/(60*60)
RIGGSlon = - RIGGSlon * RIGGS[:,9];  # RIGGS[:,9] is +1 if W, -1 if E

# throw out the ones which are not in model domain; 132 (131?) remain
cbar_masked.putmask(-20)
cRIGGS = tri.nn_interpolator(cbar_masked.flat)(RIGGSlon, RIGGSlat)
rig = RIGGS[cRIGGS > 0]
riglon = RIGGSlon[cRIGGS > 0]
riglat = RIGGSlat[cRIGGS > 0]

# add markers for RIGGS points, then quiver observed velocities
plot(riglon, riglat, '.k')
rigu = sin((pi/180)*rig[:,12]) * rig[:,10]
rigv = cos((pi/180)*rig[:,12]) * rig[:,10]
quiver(riglon, riglat, rigu, rigv, color='black')

# quiver the computed velocities at the same points; note reversal of u,v in model
uATrig = tri.nn_interpolator(vbar.flat)(riglon, riglat)
vATrig = tri.nn_interpolator(ubar.flat)(riglon, riglat)
quiver(riglon, riglat, uATrig, vATrig, color='red')
axis([-5.26168, 3.72207, -13, -5.42445])
xlabel('RIGGS grid longitude (deg E)'); ylabel('RIGGS grid latitude (deg N)')
title("""Color is speed in m/a.\n Arrows are observed (black) and computed 
(red) velocities at RIGGS points.""")
savefig("rossquiver.png")
show()

# report results comparable to Table 1 in (MacAyeal et al 1996)
ChiSqrActual = sum( ((uATrig - rigu)**2 + (vATrig - rigv)**2) / (30**2) )
print "chi^2 = %f" % (ChiSqrActual * (156.0/132.0))
print "maximum computed ice shelf speed is  %f." % (cbar.max())

# show observed versus computed scatter plot as in Figure 2 in (MacAyeal et al 1996)
figure(2);clf();hold(True)
plot(sqrt(uATrig**2 + vATrig**2), sqrt(rigu**2 + rigv**2), '.k')
plot([0, 1000],[0, 1000], 'k')
xlabel('PISM computed speed (m/a)'); ylabel('RIGGS observed speed (m/a)')
savefig("rossscatter.png")
show()

