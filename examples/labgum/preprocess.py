#!/usr/bin/env python

# Copyright (C) 2013, 2014, 2016, 2018 the PISM Authors

# This script sets up the bootstrap file.
# See also preprocess.sh.

import sys
import time
import subprocess
import argparse
import shlex
import numpy as np

try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

parser = argparse.ArgumentParser(description='Preprocess for validation using constant flux experiment from Sayag & Worster (2013).  Creates PISM-readable bootstrap file and a configuration overrides file.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-Mx', default=52,
                    help='number of points in each direction on a square grid; note MX -> cell width cases: 52 -> 10mm,  104 -> 5mm, 209 -> 2.5mm, 520 -> 1mm')
parser.add_argument('-o', metavar='FILENAME', default='initlab52.nc',
                    help='output file name to create (NetCDF)')
args = parser.parse_args()


def create_config():
    print("  creating PISM-readable config override file gumparams.nc ...")
    nc = CDF("gumparams.nc", 'w')
    config = nc.createVariable("pism_overrides", 'i4')

    attrs = {
        "constants.standard_gravity": 9.81,
        "constants.standard_gravity_doc": "m s-2; = g",

        "constants.ice.density": 1000.0,
        "constants.ice.density_doc": "kg m-3; 1% Xanthan gum in water has same density as water",

        "stress_balance.sia.bed_smoother.range": -1.0,
        "stress_balance.sia.bed_smoother.range_doc": "m; negative value de-activates bed smoother",

        "bootstrapping.defaults.geothermal_flux": 0.0,
        "bootstrapping.defaults.geothermal_flux_doc": "W m-2; no geothermal",

        "output.runtime.time_unit_name": "second",
        "output.runtime.time_unit_name_doc": "stdout uses seconds (not years) to show model time",

        "output.runtime.time_use_calendar": "no",
        "output.runtime.time_use_calendar_doc": "stdout does not use a calendar to show model time",

        "output.runtime.volume_scale_factor_log10": -15,
        "output.runtime.volume_scale_factor_log10_doc": "; an integer; log base 10 of scale factor to use for volume in summary line to stdout; -15 gives volume in cm^3",

        "output.runtime.area_scale_factor_log10": -10,
        "output.runtime.area_scale_factor_log10_doc": "; an integer; log base 10 of scale factor to use for area in summary line to stdout; -10 gives area in cm^2",

        "geometry.ice_free_thickness_standard": 1e-8,
        "geometry.ice_free_thickness_standard_doc": "m; only if the fluid is less than this is a cell marked as ice free",

        "output.ice_free_thickness_standard": 1e-8,
        "output.ice_free_thickness_standard_doc": "fluid layer exceeding this thickness is included in computations of areas and volumes",

        "time_stepping.adaptive_ratio": 0.08,
        "time_stepping.adaptive_ratio_doc": "; compare default 0.12; needs to be smaller because gum suspension is more shear-thinning than ice?",

        "stress_balance.sia.Glen_exponent": 5.9,
        "stress_balance.sia.Glen_exponent_doc": "; : n;  Sayag & Worster (2013) give n = 5.9 +- 0.2",

        "flow_law.isothermal_Glen.ice_softness": 9.7316e-09,  # vs (e.g.) 4e-25 Pa-3 s-1 for ice
        "ice_softness_doc": "Pa-n s-1; = A_0 = B_0^(-n) = (2 x 11.4 Pa s^(1/n))^(-n);  Sayag & Worster (2013) give B_0/2 = tilde mu = 11.4 +- 0.25 Pa s^(1/n)"
    }

    keys = list(attrs.keys())
    keys.sort()
    for k in keys:
        config.setncattr(k, attrs[k])

    nc.close()


create_config()

# lab setup is table with hole in the middle into which is piped the
# shear-thinning fluid, which is Xanthan gum 1% solution
Lx = 260.0e-3    # m;  = 260 mm;  maximum observed radius is 25.2 cm so we go out just a bit
Ly = Lx          # square table
flux = 3.8173e-3  # kg s-1;  = 3.8173 g s-1; Sayag personal communication
pipeR = 8.0e-3   # m;  = 8 mm;  input pipe has this radius; Sayag personal communication
temp = 20.0      # C;  fluid is at 20 deg (though it should not matter)

# set up the grid:
Mx = int(args.Mx)
My = Mx
print("  creating grid of Mx = %d by My = %d points ..." % (Mx, My))
dx = (2.0 * Lx) / float(Mx)
dy = (2.0 * Ly) / float(My)
print("  cells have dimensions dx = %.3f mm by dy = %.3f mm ..." % (dx * 1000.0, dy * 1000.0))
x = np.linspace(-Lx - dx / 2.0, Lx + dx / 2.0, Mx)
y = np.linspace(-Ly - dy / 2.0, Ly + dy / 2.0, My)

# create dummy fields
[xx, yy] = np.meshgrid(x, y)  # if there were "ndgrid" in numpy we would use it

topg = np.zeros((Mx, My))
thk = np.zeros((Mx, My))  # no fluid on table at start
artm = np.zeros((Mx, My)) + 273.15 + temp  # 20 degrees Celsius

# smb = flux as m s-1, but scaled so that the total is correct even on a coarse grid
smb = np.zeros((Mx, My))
smb[xx ** 2 + yy ** 2 <= pipeR ** 2 + 1.0e-10] = 1.0
smbpos = sum(sum(smb))
if smbpos == 0:
    print("gridding ERROR: no cells have positive input flux ... ending now")
    sys.exit(1)
else:
    print("  input flux > 0 at %d cells ..." % smbpos)
smb = (flux / (smbpos * dx * dy)) * smb  # [flux] = kg s-1  so now  [smb] = kg m-2 s-1

# Write the data:
nc = CDF(args.o, "w", format='NETCDF3_CLASSIC')  # for netCDF4 module

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

fill_value = np.nan


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

smb_var = def_var(nc, "climatic_mass_balance", "kg m-2 s-1", fill_value)
smb_var.standard_name = "land_ice_surface_specific_mass_balance_flux"
smb_var[:] = smb

artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
artm_var[:] = artm

# set global attributes
nc.Conventions = 'CF-1.4'
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()

print('  ... PISM-bootable NetCDF file %s written' % args.o)
