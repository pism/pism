#!/usr/bin/env python
#
# Copyright (C) 2011, 2014, 2018 Andy Aschwanden
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


import sys
import time
import numpy as np
from pyproj import Proj
from sys import stderr

write = stderr.write

from netCDF4 import Dataset as CDF

from optparse import OptionParser

ice_density = 910.0

__author__ = "Andy Aschwanden"

# Create PISM-readable input file from Storglaciaren DEM

parser = OptionParser()
parser.usage = "usage: %prog [options]"
parser.description = "Preprocess Storglaciaren files."


(options, args) = parser.parse_args()


# Create PISM-readable input file from Storglaciaren DEM


write('------------------------------\n')
write('PISM-Storglaciaren example\n')
write('------------------------------\n')

# data dir
data_dir = './'
# X is Northing (http://en.wikipedia.org/wiki/Swedish_grid)
XFile = data_dir + 'X.txt.gz'
# Y is Easting (http://en.wikipedia.org/wiki/Swedish_grid)
YFile = data_dir + 'Y.txt.gz'
# Bed and Surface DEMs for Storglaciaren
zBaseFile = data_dir + 'zBase.txt.gz'
zSurfFile = data_dir + 'zSurf.txt.gz'

# load coordinate information. Note: Swedish grid (RT90) uses inverse notation
# X -> northing, Y -> easting
try:
    write('Reading northing coordinate infos from %s: ' % XFile)
    X = np.loadtxt(XFile)
    write('Done.\n')
    write('Reading easting coordinate infos from %s: ' % YFile)
    Y = np.loadtxt(YFile)
    write('Done.\n')
except IOError:
    write('ERROR: File %s or %s could not be found.\n' % XFile % YFile)
    exit(2)

# load Bed DEM
try:
    write('Reading DEM from %s: ' % zBaseFile)
    zBase = np.loadtxt(zBaseFile)
    write('Done.\n')
except IOError:
    write('ERROR: File %s could not be found.\n' % zBaseFile)
    exit(2)

# load Surface DEM
try:
    write('Reading DEM from %s: ' % zSurfFile)
    zSurf = np.loadtxt(zSurfFile)
    write('Done.\n')
except IOError:
    write('ERROR: File %s could not be found.\n' % zSurfFile)
    exit(2)

# Grid size. DEM has 10m spacing.
N = zBase.shape[1] - 1
M = zBase.shape[0] - 1
Ne = N + 30
Me = M
n0 = X[-1, 0]
e0 = Y[0, 0]
de = 10  # m
dn = 10  # m

e1 = e0 + (Ne - 1) * de
n1 = n0 + (Me - 1) * dn

easting = np.linspace(e0, e1, Ne)
northing = np.linspace(n0, n1, Me)

# convert to lat/lon

ee, nn = np.meshgrid(easting, northing)

projection = Proj(init="epsg:3021")
longitude, latitude = projection(ee, nn, inverse=True)

write("Coordinates of the lower-left grid corner:\n"
      "  easting  = %.0f\n"
      "  northing = %.0f\n"
      "Grid size:\n"
      "  rows = %d\n"
      "  columns = %d\n" % (e0, n0, N, M))

# Fill value
fill_value = -9999
bed_valid_min = -5000.0
thk_valid_min = 0.0

bed = np.flipud(zBase)
dem = np.flipud(zSurf)  # ignored by bootstrapping
thk = np.flipud(zSurf-zBase)  # used for bootstrapping

m_bed = np.zeros((Me, Ne)) + 1100
m_bed[0:M, 0:N] = bed[0:M, 0:N]

m_dem = np.zeros((Me, Ne))
m_dem[0:M, 0:N] = dem[0:M, 0:N]

m_thk = np.zeros((Me, Ne))
m_thk[0:M, 0:N] = thk[0:M, 0:N]

# Replace NaNs with zeros
thk = np.nan_to_num(thk)
m_thk = np.nan_to_num(m_thk)
dem = np.nan_to_num(dem)
m_dem = np.nan_to_num(m_dem)

# There are some negative thickness values
# Quick and dirty: set to zero
# some inconsistencies in the original data still needs to be sorted out
# (filtering)
thk[thk < 0] = 0
m_thk[m_thk < 0] = 0


# Output filename
ncfile = 'pism_storglaciaren_3d.nc'

# Write the data:
nc = CDF(ncfile, "w", format='NETCDF3_CLASSIC')  # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=easting.shape[0])
nc.createDimension("y", size=northing.shape[0])

x = nc.createVariable("x", 'f4', dimensions=("x",))
x.units = "m"
x.long_name = "easting"
x.standard_name = "projection_x_coordinate"

y = nc.createVariable("y", 'f4', dimensions=("y",))
y.units = "m"
y.long_name = "northing"
y.standard_name = "projection_y_coordinate"

x[:] = easting
y[:] = northing

from scipy.interpolate import griddata


def def_var(nc, name, units, fillvalue):
    var = nc.createVariable(name, 'f', dimensions=("y", "x"), fill_value=fillvalue)
    var.units = units
    return var


lon_var = def_var(nc, "lon", "degrees_east", None)
lon_var.standard_name = "longitude"
lon_var[:] = longitude

lat_var = def_var(nc, "lat", "degrees_north", None)
lat_var.standard_name = "latitude"
lat_var[:] = latitude

bed_var = def_var(nc, "topg", "m", fill_value)
bed_var.valid_min = bed_valid_min
bed_var.standard_name = "bedrock_altitude"
bed_var.coordinates = "lat lon"
bed_var[:] = m_bed

thk_var = def_var(nc, "thk", "m", fill_value)
thk_var.valid_min = thk_valid_min
thk_var.standard_name = "land_ice_thickness"
thk_var.coordinates = "lat lon"
thk_var[:] = m_thk

dem_var = def_var(nc, "usurf_from_dem", "m", fill_value)
dem_var.standard_name = "surface_altitude"
dem_var.coordinates = "lat lon"
dem_var[:] = m_dem

ftt_mask = def_var(nc, "ftt_mask", "", fill_value)
ftt_mask[:] = 1


# generate (somewhat) reasonable acab
acab_max = 2.5  # m/a
acab_min = -3.0  # m/a
acab_up = easting.min() + 200  # m; location of upstream end of linear acab
acab_down = easting.max() - 900  # m;location of downstream end of linear acab
acab = np.ones_like(longitude)

acab[:] = acab_min
acab[:] = acab_max - (acab_max-acab_min) * (easting - acab_up) / (acab_down - acab_up)
acab[m_thk < 1] = acab_min

acab_var = def_var(nc, "climatic_mass_balance", "kg m-2 year-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = acab * ice_density

# Set boundary conditions for Scandinavian-type polythermal glacier
# ------------------------------------------------------------------------------
#
# (A) Surface temperature for temperature equation bc
T0 = 273.15  # K
Tma = -6.0  # degC, mean annual air temperature at Tarfala
zcts = 1300   # m a.s.l.; altitude where CTS is at the surface, projected to topg
slope = 100    # m; range around which surface temp transition happens


# smoothed version; FIXME:  can't we at least have it depend on initial DEM?
#   additional lapse rate?
artm = T0 + Tma * (zcts + slope - m_bed) / (2.0 * slope)
artm[m_bed < zcts-slope] = T0 + Tma
artm[m_bed > zcts+slope] = T0

artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
artm_var[:] = artm

# set global attributes
nc.Conventions = "CF-1.4"
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)
# proj4 string equivalent to EPSG 3021
nc.projection = "+proj=tmerc +lat_0=0 +lon_0=15.80827777777778 +k=1 +x_0=1500000 +y_0=0 +ellps=bessel +towgs84=419.384,99.3335,591.345,0.850389,1.81728,-7.86224,-0.99496 +units=m +no_defs"
nc.close()
write('Done writing NetCDF file %s!\n' % ncfile)
