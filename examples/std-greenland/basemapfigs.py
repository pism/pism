#!/usr/bin/env python

# generate figures in Getting Started section of User's Manual

# usage:
#   $ python basemapfigs.py FILEROOT [FIELD] [DPI]
# where
#   FILEROOT   root of NetCDF filename and output .png figures
#   FIELD      optional: one of {velbase_mag, [velsurf_mag], mask, usurf}  (edit script to add more)
#   DPI        optional: resolution in dots per inch [200]
#
# equivalent usages:
#   $ python basemapfigs.py g20km_10ka_hy velsurf_mag 200
#   $ python basemapfigs.py g20km_10ka_hy velsurf_mag
#   $ python basemapfigs.py g20km_10ka_hy
#
# generate figs like those in Getting Started section of User's Manual:
#   $ for FLD in velsurf_mag usurf velbase_mag mask; do ./basemapfigs.py g20km_10ka_hy ${FLD}; done
#
# crop out western Greenland with command like this (uses ImageMagick):
#   $ ./basemapfigs.py g20km_10ka_hy velsurf_mag 500
#   $ convert -crop 600x800+400+800 +repage g20km_10ka_hy-velsurf_mag.png g20km-detail.png
#
# batch generate figures from a parameter study like this:
#   $ for QQ in 0.1 0.5 1.0; do for EE in 1 3 6; do ../basemapfigs.py p10km_q${QQ}_e${EE} velsurf_mag 100; done; done
#   $ for QQ in 0.1 0.5 1.0; do for EE in 1 3 6; do convert -crop 274x486+50+6 +repage p10km_q${QQ}_e${EE}-velsurf_mag.png p10km-${QQ}-${EE}-csurf.png; done; done

from mpl_toolkits.basemap import Basemap

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import sys

if len(sys.argv) < 2:
    print("ERROR: first argument must be root of filename ...")
    sys.exit(1)
rootname = sys.argv[1]
try:
    nc = NC(rootname + '.nc', 'r')
except:
    print("ERROR: can't read from file %s.nc ..." % rootname)
    sys.exit(2)

if len(sys.argv) >= 3:
    field = sys.argv[2]
else:
    field = 'velsurf_mag'

if len(sys.argv) >= 4:
    mydpi = float(sys.argv[3])
else:
    mydpi = 200

bluemarble = False  # if True, use Blue Marble background

if (field == 'velsurf_mag') | (field == 'velbase_mag'):
    fill = nc.variables[field]._FillValue
    logscale = True
    contour100 = True
    myvmin = 1.0
    myvmax = 6.0e3
    ticklist = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
elif field == 'surfvelmag':
    fill = 0.0
    logscale = True
    contour100 = True
    myvmin = 1.0
    myvmax = 6.0e3
    ticklist = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
elif field == 'usurf':
    fill = 0.0
    logscale = False
    contour100 = False
    myvmin = 1.0
    myvmax = 3500.0
    ticklist = [100, 500, 1000, 1500, 2000, 2500, 3000, 3500]
elif field == 'mask':
    fill = -1.0
    logscale = False
    contour100 = False
    myvmin = 0.0
    myvmax = 4.0
    ticklist = [0, 1, 2, 3, 4]
elif field == 'basal_melt_rate_grounded':
    fill = -2.0e+09
    logscale = True
    contour100 = False
    myvmin = 0.9e-4
    myvmax = 1.1
    ticklist = [0.0001, 0.001, 0.01, 0.1, 1.0]
elif field == 'tillwat':
    fill = -2.0e+09
    logscale = False
    contour100 = False
    myvmin = 0.0
    myvmax = 2.0
    ticklist = [0.0, 0.5, 1.0, 1.5, 2.0]
elif field == 'bwat':
    fill = -2.0e+09
    logscale = True
    contour100 = False
    myvmin = 0.9e-4
    myvmax = 1.1
    ticklist = [0.0001, 0.001, 0.01, 0.1, 1.0]
elif field == 'bwprel':
    fill = -2.0e+09
    logscale = False
    contour100 = False
    myvmin = 0.0
    myvmax = 1.0
    ticklist = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
else:
    print('invalid choice for FIELD option')
    sys.exit(3)

# we need to know longitudes and latitudes corresponding to grid
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
if field == 'surfvelmag':
    lon = np.squeeze(lon).transpose()
    lat = np.squeeze(lat).transpose()

# x and y *in the dataset* are only used to determine plotting domain
# dimensions
if field == 'surfvelmag':
    x = nc.variables['x1'][:]
    y = nc.variables['y1'][:]
else:
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
width = x.max() - x.min()
height = y.max() - y.min()

# load data
if field == 'bwprel':
    thkvar = np.squeeze(nc.variables['thk'][:])
    myvar = np.squeeze(nc.variables['bwp'][:])
    myvar = np.ma.array(myvar, mask=(thkvar == 0.0))
    thkvar = np.ma.array(thkvar, mask=(thkvar == 0.0))
    myvar = myvar / (910.0 * 9.81 * thkvar)
else:
    myvar = np.squeeze(nc.variables[field][:])

# mask out ice free etc.; note 'mask' does not get masked
if (field == 'surfvelmag'):
    myvar = myvar.transpose()
    thkvar = np.squeeze(nc.variables['thk'][:]).transpose()
    myvar = np.ma.array(myvar, mask=(thkvar == 0.0))
elif (field != 'mask'):
    maskvar = np.squeeze(nc.variables['mask'][:])
    if (field == 'basal_melt_rate_grounded') | (field == 'bwat'):
        myvar[myvar < myvmin] = myvmin
    if (field == 'usurf'):
        myvar = np.ma.array(myvar, mask=(maskvar == 4))
    else:
        myvar = np.ma.array(myvar, mask=(maskvar != 2))

m = Basemap(width=1.1 * width,    # width in projection coordinates, in meters
            height=1.05 * height,      # height
            resolution='l',     # coastline resolution, can be 'l' (low), 'h'
                                # (high) and 'f' (full)
            projection='stere',  # stereographic projection
            lat_ts=71,          # latitude of true scale
            lon_0=-41,          # longitude of the plotting domain center
            lat_0=72)           # latitude of the plotting domain center

# m.drawcoastlines()

# draw the Blue Marble background (requires PIL, the Python Imaging Library)
if bluemarble:  # seems to reverse N and S
    m.bluemarble()

# convert longitudes and latitudes to x and y:
xx, yy = m(lon, lat)

if contour100:
    # mark 100 m/a contour in black:
    m.contour(xx, yy, myvar, [100], colors="black")

# plot log color scale or not
if logscale:
    m.pcolormesh(xx, yy, myvar,
                 norm=colors.LogNorm(vmin=myvmin, vmax=myvmax))
else:
    m.pcolormesh(xx, yy, myvar, vmin=myvmin, vmax=myvmax)

# add a colorbar:
plt.colorbar(extend='both',
             ticks=ticklist,
             format="%d")

# draw parallels and meridians
# labels kwarg is where to draw ticks: [left, right, top, bottom]
m.drawparallels(np.arange(-55., 90., 5.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-120., 30., 10.), labels=[0, 0, 0, 1])

outname = rootname + '-' + field + '.png'
print("saving image to file %s ..." % outname)
plt.savefig(outname, dpi=mydpi, bbox_inches='tight')
