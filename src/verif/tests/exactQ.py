#!/usr/bin/env python

# Computes and plots exact similarity solution "test Q".  See equation
# (3.19) in [\ref PeglerListerWorster2012] = PLW2012:
#   S. Pegler, J. Lister, and M. G. Worster, 2012.  Release of a viscous
#   power-law fluid over an inviscid ocean", J. Fluid Mech. 700, 63--76.

from pylab import *
import sys
import time
# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

#   $ ./exactQ.py -o foo.png ts_lab51.nc

from optparse import OptionParser

parser = OptionParser()
parser.usage = "usage: %prog Mx OUTFILE"
parser.description = "A script which builds bootstrap file for Test Q."

(options, args) = parser.parse_args()

nts = len(args)
if (nts<2) | (nts>2):
  print "ERROR; needs two arguments:  exactQ.py Mx OUTFILE"
  print "  ... EXITING"
  exit(-1)

Mx = int(args[0])
My = Mx

SperA= 31556926.0

g    = 9.81
rho  = 910.0    # density of ice; kg/m^3
rhow = 1028.0   # density of ocean water; kg/m^3

n = 3.0
barB = 1.9e8    # strength of shelf; Pa s^(1/3); from MacAyeal et al 1996;
                # is this equal to \tilde mu or not?

def timeQ(H0,R0):
  m = (1.0/n) - 1.0
  gprime = (rhow - rho) * g / rhow     # see just after (2.1) in PLW2012
  nurescale = 3.0**(m/2) * barB / rho  # see just after (3.19) in PLW2012
  C0 = 12.0 * nurescale / gprime       # (3.19)
  t = (1.0 / (2.0 * n)) * (C0 / H0)**n # (3.19)
  return t

# similarity solution: choose dimensions, get told the time and volume
H0 = 1000.0     # m
R0 = 100.0e3    # m; 100 km
t = timeQ(H0,R0)
V = pi * R0**2 * H0

print 'exact Test Q has the following parameters for the start time:'
print '  time      t  = %.3e s = %f a' % (t, t/SperA)
print '  volume    V  = %.3e m^3 = %.3e km^3' % (V, V/1.0e9)
print '  thickness H0 = %.3f m' % H0
print '  radius    R0 = %.3f m = %.3f km' % (R0, R0/1.0e3)

ncfile = args[1]

print 'building PISM-bootstrapable file %s with Mx = %d and My = %d grid points ...' % (ncfile,Mx,My)

# set up the grid:
Lx = 200.0e3
Ly = Lx
x = linspace(-Lx,Lx,Mx)
y = linspace(-Ly,Ly,My)
[xx,yy] = meshgrid(x,y)         # if there were "ndgrid" in numpy?
rr = sqrt(xx**2 + yy**2)

# create fields
#floatdepth = (rho / rhow) * H0
#topg = zeros((Mx,My)) - depth + 10.0 - ...  # so it is floating
topg = zeros((Mx,My)) - 1200.0  # so it is floating
artm = zeros((Mx,My)) + 263.15  # -10 degrees Celsius; just a guess
acab = zeros((Mx,My))
thk  = zeros((Mx,My))
thk[xx**2 + yy**2 <= R0**2] = H0;
zerossabc = zeros((Mx,My))

# Write the data:
nc = CDF(ncfile, "w",format='NETCDF3_CLASSIC') # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=Mx)
nc.createDimension("y", size=My)

x_var = nc.createVariable("x", 'f4', dimensions=("x",))
x_var.units = "m";
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = x

y_var = nc.createVariable("y", 'f4', dimensions=("y",))
y_var.units = "m";
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"
y_var[:] = y

fill_value = nan

def def_var(nc, name, units, fillvalue):
    # dimension transpose is standard: "float thk(y, x)" in NetCDF file
    var = nc.createVariable(name, 'f', dimensions=("y", "x"), fill_value=fillvalue)
    var.units = units
    return var

bed_var = def_var(nc, "topg", "m", fill_value)
bed_var.standard_name = "bedrock_altitude"
bed_var[:] = topg

thk_var = def_var(nc, "thk", "m", fill_value)
thk_var.standard_name = "land_ice_thickness"
thk_var[:] = thk

acab_var = def_var(nc, "climatic_mass_balance", "m s-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = acab

artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
artm_var[:] = artm

u_ssa_bc_var = def_var(nc, "u_ssa_bc", "m s-1", fill_value)
u_ssa_bc_var[:] = zerossabc.copy()

v_ssa_bc_var = def_var(nc, "v_ssa_bc", "m s-1", fill_value)
v_ssa_bc_var[:] = zerossabc.copy()

bcflag_var = nc.createVariable("bcflag", "i", dimensions=("y", "x"))
bcflag_var[:] = ((xx == 0.0) & (yy == 0.0))


# set global attributes
nc.Conventions = 'CF-1.4'
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()
print('PISM-bootstrappable NetCDF file %s written' % ncfile)
print('  ... now run   FIXME')


# don't use  -pik  because it will do  -kill_icebergs  and the whole thing is an iceberg

# could use  -ssa_view_nuh  if desired

# trial run:
#   $ MX=101
#   $ ./exactQ.py $MX initQ$MX.nc
#   $ pismr -o outQ$MX.nc -boot_file initQ$MX.nc -Mx $MX -My $MX -Mz 21 -Lz 1500 -z_spacing equal -surface given -no_sia -no_energy -ssa_floating_only -ssa_dirichlet_bc -cfbc -part_grid -part_redist -y 0 -o_order zyx -ssa_e 1.0 -ssa_flow_law isothermal_glen

