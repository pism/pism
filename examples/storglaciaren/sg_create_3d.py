#!/usr/bin/env python
#
# Copyright (C) 2011 Andy Aschwanden
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
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

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

from optparse import OptionParser

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
# Bed and Surface DEMs for Storglaciaren
XFile     = data_dir + 'X.txt.gz'
YFile     = data_dir + 'Y.txt.gz'
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
    write('ERROR: File %s or %s could not be found.\n'  % XFile % YFile)
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
N  = zBase.shape[1]
M  = zBase.shape[0]
e0 = Y.min()
n0 = X.min()
de = 10 # m
dn = 10 # m

e1 = e0 + (N-1)*de
n1 = n0 + (M-1)*dn

easting  = np.linspace(e0, e1, N)
northing = np.linspace(n0, n1, M)

# convert to lat/lon

# From http://lists.maptools.org/pipermail/proj/2008-December/004165.html:
#
# However, a simpler method, now recommended by the Swedish Land Survey
# instead of a 7-parameter shift, is to start from the WGS84 datum, and than
# tweak the projection parameters a little: just use a Transverse Mercator
# with
#   central meridian: 15" 48' 22.624306" E  
#   scale factor:     1.00000561024
#   false easting:    1500064.274 m
#   false northing:   -667.711 m
# ( http://www.lantmateriet.se/templates/LMV_Page.aspx?id=5197&lang=EN )

projRT90 = "+proj=tmerc +datum=WGS84 +lon_0=-15.806284 +x_0=1500064.274 +y_0=-667.711 +k=1.00000561024 +units=m"

ee, nn = np.meshgrid(easting, northing)

projection = Proj(projRT90)
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
dem = np.flipud(zSurf) # ignored by bootstrapping
thk = np.flipud(zSurf-zBase) # used for bootstrapping

# Replace NaNs with zeros
thk = np.nan_to_num(thk)
# There are some negative thickness values
# Quick and dirty: set to zero
# some inconsistencies in the original data still needs to be sorted out
# (filtering)
thk[thk<0] = 0



# Output filename
ncfile = 'pism_storglaciaren_3d.nc'

# Write the data:
nc = CDF(ncfile, "w",format='NETCDF3_CLASSIC') # for netCDF4 module

# Create dimensions x and y
nc.createDimension("x", size=easting.shape[0])
nc.createDimension("y", size=northing.shape[0])

x = nc.createVariable("x", 'f4', dimensions=("x",))
x.units = "m";
x.long_name = "easting"
x.standard_name = "projection_x_coordinate"

y = nc.createVariable("y", 'f4', dimensions=("y",))
y.units = "m";
y.long_name = "northing"
y.standard_name = "projection_y_coordinate"

x[:] = easting
y[:] = northing


def def_var(nc, name, units, fillvalue):
    var = nc.createVariable(name, 'f', dimensions=("y", "x"),fill_value=fillvalue)
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
bed_var[:] = bed

thk_var = def_var(nc, "thk", "m", fill_value)
thk_var.valid_min = thk_valid_min
thk_var.standard_name = "land_ice_thickness"
thk_var.coordinates = "lat lon"
thk_var[:] = thk

dem_var = def_var(nc, "usurf_from_dem", "m", fill_value)
dem_var.standard_name = "surface_altitude"
dem_var.coordinates = "lat lon"
dem_var[:] = dem

# generate (somewhat) reasonable acab
acab_max =  2.5 # m/a
acab_min = -3.0 # m/a
acab_up = easting.min() + 200 # m; location of upstream end of linear acab
acab_down = easting.max() - 600 # m;location of downstream end of linear acab
acab = np.ones_like(dem)

acab[:] = acab_max - (acab_max-acab_min) * (easting - acab_up) / (acab_down - acab_up)
acab[thk<1] = acab_min

acab_var = def_var(nc, "climatic_mass_balance", "m year-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = acab

# Set boundary conditions for Scandinavian-type polythermal glacier
# ------------------------------------------------------------------------------
#
# (A) Surface temperature for temperature equation bc
T0    = 273.15 # K
Tma   =  -6.0  # degC, mean annual air temperature at Tarfala
zcts  = 1300   # m a.s.l.; altitude where CTS is at the surface, projected to topg
slope  = 100    # m; range around which surface temp transition happens

# old abrupt jump:
#artm  = np.zeros((M,N),float) + T0
#artm[bed<zcts] = T0 + Tma # Scandinavian-type polythermal glacier

# smoothed version; FIXME:  can't we at least have it depend on initial DEM?
#   additional lapse rate?
artm = T0 + Tma * (zcts + slope - bed) / (2.0 * slope)
artm[bed<zcts-slope] = T0 + Tma
artm[bed>zcts+slope] = T0

artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
artm_var[:] = artm

# set global attributes
nc.Conventions = "CF-1.4"
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)
nc.projection = projRT90
nc.close()
write('Done writing NetCDF file %s!\n' % ncfile)

