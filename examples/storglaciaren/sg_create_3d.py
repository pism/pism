#!/usr/bin/env python
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

# default values
ISTEMP=False

parser = OptionParser()
parser.usage = "usage: %prog [options]"
parser.description = "Preprocess Storglaciaren files."
parser.add_option("-t", "--temperate",dest="istemp",action="store_true",
                  help="sets upper surface boundary condition to temperate",default=ISTEMP)


(options, args) = parser.parse_args()
istemp = options.istemp


# Create PISM-readable input file from Storglaciaren DEM
# 05/13/09 AA started working on script
# 05/14/09 AA reads and writes fields
# 05/15/09 AA input field moved to ../data/
# 05/20/09 AA x and y uses the Swedish grid (or should this go to lat/lon?)
# 06/09/09 AA corrected handling of x,y coordinates, acab variable added


write('UAF-Storglaciaren PISM project\n')
write('------------------------------\n')

# data dir
data_dir = '../data/'
# Bed and Surface DEMs for Storglaciaren
XFile     = data_dir + 'X.txt'
YFile     = data_dir + 'Y.txt'
zBaseFile = data_dir + 'zBase.txt'
zSurfFile = data_dir + 'zSurf.txt'
MBFile    = data_dir + 'mb2001.txt'

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

# load Annual Mass Balance 2001 (in meters)
try:
    write('Reading annual mass balance (2001) from %s: ' % MBFile)
    MB = np.loadtxt(MBFile)
    write('Done.\n')
except IOError:
    write('ERROR: File %s could not be found.\n' % MBFile)
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
      "  easting  = %f\n"
      "  northing = %f\n"
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
mb  = MB # no flip required

# Replace NaNs with zeros
thk = np.nan_to_num(thk)
# There are some negative thickness values
# Quick and dirty: set to zero
# some inconsistencies in the original data still needs to be sorted out
# (filtering)
thk[thk<0] = 0

# insert mild ablation where no data are available, off the glacier
mb[thk<1] = -3.

# Set boundary conditions
# ------------------------------------------------------------------------------
#
# (A) Surface temperature for temperature equation bc
T0    = 273.15 # K
Tma   =  -6.0  # degC, mean annual air temperature at Tarfala
zcts  = 1300   # m a.s.l.; altitude where CTS is at the surface, projected to topg
zbts  = 1250   # m a.s.l.; altitude where CTS is at the base; just a wild guess

if istemp:
    artm  = np.zeros((M,N),float) + T0
else:
    artm  = np.zeros((M,N),float) + T0
    artm[bed<zcts] = T0 + Tma # Scandinavian-type polythermal glacier


# Output filename
ncfile = 'pism_storglaciaren.nc'

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
    var = nc.createVariable(name, 'f', dimensions=("y", "x"))
    var.units = units
    if (fillvalue != None):
        var._FillValue = fillvalue
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
bed_var[:] = bed

thk_var = def_var(nc, "thk", "m", fill_value)
thk_var.valid_min = thk_valid_min
thk_var.standard_name = "land_ice_thickness"
thk_var[:] = thk

dem_var = def_var(nc, "usurf_from_dem", "m", fill_value)
dem_var[:] = dem

acab_var = def_var(nc, "acab", "m year-1", fill_value)
acab_var.standard_name = "land_ice_surface_specific_mass_balance"
acab_var[:] = mb

artm_var = def_var(nc, "artm", "K", fill_value)
artm_var[:] = artm

# set global attributes
nc.Conventions = "CF-1.4"
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(nc, 'history', historystr)

nc.close()
write('Done writing NetCDF file %s!\n' % ncfile)

np.savetxt('acab2001.txt',mb)
np.savetxt('thk.txt',thk)
np.savetxt('topg.txt',bed)
        
