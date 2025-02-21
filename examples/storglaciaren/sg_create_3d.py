#!/usr/bin/env python3
#
# Copyright (C) 2025 Andy Aschwanden
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
import xarray as xr
import rioxarray

ice_density = 910.0

__author__ = "Andy Aschwanden"

# Create PISM-readable input file from Storglaciaren DEM


print('------------------------------\n')
print('PISM-Storglaciaren 3D example\n')
print('------------------------------\n')

# data dir
data_dir = './'
# X is Northing (http://en.wikipedia.org/wiki/Swedish_grid)
XFile = data_dir + 'X.txt.gz'
# Y is Easting (http://en.wikipedia.org/wiki/Swedish_grid)
YFile = data_dir + 'Y.txt.gz'
# Bed and Surface DEMs for Storglaciaren
zBaseFile = data_dir + 'zBase.txt.gz'
zSurfFile = data_dir + 'zSurf.txt.gz'

print('Reading northing coordinate infos from %s: ' % XFile)
X = np.loadtxt(XFile)
print('Done.\n')
print('Reading easting coordinate infos from %s: ' % YFile)
Y = np.loadtxt(YFile)
print('Done.\n')

# load Bed DEM
print('Reading DEM from %s: ' % zBaseFile)
zBase = np.loadtxt(zBaseFile)
print('Done.\n')


print('Reading DEM from %s: ' % zSurfFile)
zSurf = np.loadtxt(zSurfFile)
print('Done.\n')

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

ee, nn = np.meshgrid(easting, northing)


print("Coordinates of the lower-left grid corner:\n"
      "  easting  = %.0f\n"
      "  northing = %.0f\n"
      "Grid size:\n"
      "  rows = %d\n"
      "  columns = %d\n" % (e0, n0, N, M))

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
m_dem[m_thk==0] = m_bed[m_thk==0]

x: xr.DataArray = xr.DataArray(easting, coords={"x": easting.shape[0]}, attrs={"units": "m", "_FillValue": False}, name="x")
y: xr.DataArray = xr.DataArray(northing, coords={"y": northing.shape[0]}, attrs={"units": "m", "_FillValue":False}, name="y")

bed = xr.DataArray(m_bed, dims=["y", "x"], coords={"y": y, "x": x}, attrs={"standard_name": "bedrock_altitude", "units": "m"}, name="bed")
thickness = xr.DataArray(m_thk, dims=["y", "x"], coords={"y": y, "x": x}, attrs={"standard_name": "land_ice_thickness", "units": "m"}, name="thickness")
usurf = xr.DataArray(m_dem, dims=["y", "x"], coords={"y": y, "x": x}, attrs={"standard_name": "surface_altitude", "units": "m"}, name="surface")

ftt_mask = 0 * (m_thk>0) + 1 * (m_thk==0)
mask = xr.DataArray(ftt_mask, dims=["y", "x"], coords={"y": y, "x": x}, name="ftt_mask")


# generate (somewhat) reasonable acab
acab_max = 2.5  # m/a
acab_min = -3.0  # m/a
acab_up = easting.min() + 200  # m; location of upstream end of linear acab
acab_down = easting.max() - 900  # m;location of downstream end of linear acab
acab = np.ones_like(bed)

acab[:] = acab_min
acab[:] = acab_max - (acab_max-acab_min) * (easting - acab_up) / (acab_down - acab_up)
acab[m_thk < 1] = acab_min

cmb = xr.DataArray(acab * ice_density, dims=["y", "x"], coords={"y": y, "x": x},
                   attrs={"standard_name": "land_ice_surface_specific_mass_balance", "units": "kg m^-2 yr^-1"},
                   name="climatic_mass_balance")

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

temp = xr.DataArray(artm, dims=["y", "x"], coords={"y": y, "x": x}, attrs={"unit": "kelvin"}, name="ice_surface_temp")

ds = xr.merge([bed, thickness, usurf, mask, cmb, temp])
ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y")
ds.rio.write_crs("EPSG:3021", inplace=True)
ds.to_netcdf("pism_storglaciaren_3d.nc")
