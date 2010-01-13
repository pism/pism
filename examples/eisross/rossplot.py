#!/usr/bin/env python
# ross_plot.py plots the computed Ross ice shelf speed compared to observed
# values from RIGGS data.
#
# This script depends on the following Python packages:
#   1) numpy (see http://numpy.scipy.org/; python-numpy is Debian package)
#   2) matplotlib (see http://matplotlib.sourceforge.net/; python-matplotlib 
#        is Debian package)
#   3) netcdf-python (see http://code.google.com/p/netcdf4-python/; no Debian package)
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
# CK 27may08, ..., 12jan10

from numpy import ma, loadtxt, squeeze, linspace, tile, repeat, sin, pi, cos, sqrt
from pylab import figure, clf, hold, pcolor, colorbar, plot, quiver, axis, xlabel, ylabel, savefig, show
from getopt import getopt, GetoptError
from sys import argv, exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

try:
    from scikits.delaunay import Triangulation
except:
    print "ERROR:  scikits.delaunay  not installed (?)"
    print "see PISM Installation Manual for information"
    exit(1)

seconds_per_year = 3.1556926e7

# process command line arguments
try:
    opts, args = getopt(argv[1:], "p:r:", ["pism-output=", "riggs="])
    # defaults:
    pism_output = "rossComputed.nc"
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
    H = squeeze(infile.variables["thk"][:])
    mask = squeeze(infile.variables["mask"][:])
    cbar = squeeze(infile.variables["cbar"][:])
    ubar = squeeze(infile.variables["uvelsurf"][:])
    vbar = squeeze(infile.variables["vvelsurf"][:])
    print "done."
except Exception:
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

cbar_masked = cbar_masked.filled(-20)

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
axis([-5.26168, 3.72207, -13, -5.42445])
xlabel('RIGGS grid longitude (deg E)'); ylabel('RIGGS grid latitude (deg N)')
#title("""Color is speed in m/a.\n Arrows are observed (black) and computed 
#(red) velocities at RIGGS points.""")
print "saving figure ross_quiver.png"
savefig("ross_quiver.png")

# to report results comparable to Table 1 in (MacAyeal et al 1996)
#ChiSqrActual = sum( ((uATrig - rigu)**2 + (vATrig - rigv)**2) / (30**2) )
#print "chi^2 = %f" % (ChiSqrActual * (156.0/132.0))
#print "maximum computed ice shelf speed is  %f." % (cbar.max())

# show observed versus computed scatter plot as in Figure 2 in (MacAyeal et al 1996)
figure(2);clf();hold(True)
plot(sqrt(uATrig**2 + vATrig**2), sqrt(rigu**2 + rigv**2), '.k')
plot([0, 1000],[0, 1000], 'k')
xlabel('PISM computed speed (m/a)'); ylabel('RIGGS observed speed (m/a)')
print "saving figure ross_scatter.png"
savefig("ross_scatter.png")

print "pausing to show figures ..."
show()

