#!/usr/bin/env python

# Computes exact similarity solution "test Q" at an initial time and after
# a specified runtime.  Suggests a PISM run for this runtime.  The exact
# solution for velocity and for thickness can be compared to PISM output.
#
# For the similarity solution see equation (3.19) in
# [\ref PeglerListerWorster2012] = PLW2012:
#   S. Pegler, J. Lister, and M. G. Worster, 2012.  Release of a viscous
#   power-law fluid over an inviscid ocean", J. Fluid Mech. 700, 63--76.

# notes:
#   don't use  -pik  because it will do  -kill_icebergs  and the whole thing is an iceberg
#   could use  -ssa_view_nuh  if desired

# FIXME:  need to carefully set hardness

from pylab import *
import sys
import time
# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

from optparse import OptionParser

parser = OptionParser()
parser.usage = \
    """%prog MX DURATION

where MX       is number of grid points,
      DURATION is time in years for run,

Example:  Try this diagnostic-only run:
  $ export MX=101 YY=0
  $ ./exactQ.py $MX $YY
  $ pismr -o outQ$MX.nc -y $YY -i initQ$MX.nc -bootstrap -Mx $MX -My $MX -Mz 21 -Lz 1500 -z_spacing equal -surface given -stress_balance ssa -energy none -yield_stress constant -tauc 1e6 -ssa_dirichlet_bc -cfbc -part_grid -o_order zyx -ssa_e 1.0 -ssa_flow_law isothermal_glen
"""
parser.description = "A script which runs Test Q."
(options, args) = parser.parse_args()
if (len(args) < 2) | (len(args) > 2):
    print("ERROR; exactQ.py needs two arguments; run with --help to see usage")
    print("... EXITING")
    exit(-1)

SperA = 31556926.0

Mx = int(args[0])
My = Mx
runtime = float(args[1]) * SperA

ncfile = "initQ%d.nc" % Mx

# basic parameters
g = 9.81
rho = 910.0    # density of ice; kg/m^3
rhow = 1028.0   # density of ocean water; kg/m^3
n = 3.0
barB = 1.9e8    # strength of shelf; Pa s^(1/3); from MacAyeal et al 1996;
# FIXME is this equal to \tilde mu or not?

# derived parameters
m = (1.0 / n) - 1.0
gprime = (rhow - rho) * g / rhow       # see just after (2.1) in PLW2012
nurescale = 3.0 ** (m / 2.0) * barB / rho  # see just after (3.19) in PLW2012
C0 = 12.0 * nurescale / gprime         # constant in (3.19) in PLW2012


def timeQ(H0):
    # invert the first formula in (3.19) in PLW2012 to get t = f(H0)
    t = (1.0 / (2.0 * n)) * (C0 / H0) ** n
    return t


def geomvelQ(t, r, V):
    # computes (3.19) in PLW2012, assuming t,V are scalar and r is array
    # returns H,u with same size as r, but R is scalar
    tnt = 2.0 * n * t
    H = C0 * tnt ** (-1.0 / n) * ones(shape(r))
    u = r / tnt
    R = (V / (C0 * pi)) ** 0.5 * tnt ** (1.0 / (2.0 * n))
    return H, u, R


# similarity solution: choose dimensions, get told the time and volume
H0 = 1000.0     # m
R0 = 100.0e3    # m; 100 km
V = pi * R0 ** 2 * H0
t0 = timeQ(H0)

t = t0 + runtime

print('exact Test Q has the following parameters for the start time t=t0:')
print('  time      t0 = %.3e s = %f a' % (t0, t0 / SperA))
print('  thickness H0 = %.3f m' % H0)
print('  radius    R0 = %.3f km' % (R0 / 1.0e3))
print('building PISM bootstrap file %s with Mx = %d and My = %d grid points ...' % (ncfile, Mx, My))

# set up the grid:
Lx = 200.0e3
Ly = Lx
x = linspace(-Lx, Lx, Mx)
y = linspace(-Ly, Ly, My)
[xx, yy] = meshgrid(x, y)         # if there were "ndgrid" in numpy?
rr = sqrt(xx ** 2 + yy ** 2)

fill_value = nan

# create initial fields
topg = zeros((Mx, My)) - 1200.0  # so it is floating
ice_surface_temp = zeros((Mx, My)) + 263.15  # -10 degrees Celsius; just a guess
climatic_mass_balance = zeros((Mx, My))
thk = zeros((Mx, My))
thk[rr <= R0] = H0
zerossabc = zeros((Mx, My))

# create exact solution fields
thk_exact, c_exact, R_exact = geomvelQ(t, rr, V)
thk_exact[rr > R_exact] = 0.0
c_exact *= SperA
c_exact[rr > R_exact] = 0.0

print('exact Test Q at time t=%f years is in these variables:' % (t / SperA))
print('  c_exact, with max = %.3e' % c_exact.max())
print('  thk_exact, with max = %.3e' % thk_exact.max())
print('and R_exact = %.3f km' % (R_exact / 1.0e3))

# Write the data:
nc = CDF(ncfile, "w", format='NETCDF3_CLASSIC')  # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=Mx)
nc.createDimension("y", size=My)

x_var = nc.createVariable("x", 'f4', dimensions=("x",))
x_var.units = "m"
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = x

y_var = nc.createVariable("y", 'f4', dimensions=("y",))
y_var.units = "m"
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"
y_var[:] = y


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

climatic_mass_balance_var = def_var(nc, "climatic_mass_balance", "kg m-2 s-1", fill_value)
climatic_mass_balance_var.standard_name = "land_ice_surface_specific_mass_balance"
climatic_mass_balance_var[:] = climatic_mass_balance

ice_surface_temp_var = def_var(nc, "ice_surface_temp", "K", fill_value)
ice_surface_temp_var[:] = ice_surface_temp

u_ssa_bc_var = def_var(nc, "u_ssa_bc", "m s-1", fill_value)
u_ssa_bc_var[:] = zerossabc.copy()

v_ssa_bc_var = def_var(nc, "v_ssa_bc", "m s-1", fill_value)
v_ssa_bc_var[:] = zerossabc.copy()

bc_mask_var = nc.createVariable("bc_mask", "i", dimensions=("y", "x"))
bc_mask_var[:] = ((xx == 0.0) & (yy == 0.0))

thk_exact_var = def_var(nc, "thk_exact", "m", fill_value)
thk_exact_var[:] = thk_exact

c_exact_var = def_var(nc, "c_exact", "m year-1", fill_value)
c_exact_var[:] = c_exact

# set global attributes
nc.Conventions = 'CF-1.4'
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()
print('file %s written ...' % ncfile)
print('  ... now run   FIXME')
