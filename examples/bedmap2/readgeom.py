#!/usr/bin/env python

# readgeom   reads bed elevation and ice thickness from BEDMAP2 files
#            and either shows as figures or writes into NetCDF file

# edit these:
BM2PATH = '/home/ed/Desktop/bedmap2_bin/'
showgeom = False

import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt
import sys

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

from PISMNC import PISMDataset as PNC


N = 6667

print("reading bedmap2 binary files from %s ..." % (BM2PATH))

fname = BM2PATH + 'bedmap2_bed.flt'
bed = ma.masked_equal(np.reshape(np.fromfile(fname, dtype=np.float32), (N, N)), -9999.0)
fname = BM2PATH + 'bedmap2_thickness.flt'
thk = ma.masked_equal(np.reshape(np.fromfile(fname, dtype=np.float32), (N, N)), -9999.0)
fname = BM2PATH + 'bedmap2_icemask_grounded_and_shelves.flt'
msk = ma.masked_equal(np.reshape(np.fromfile(fname, dtype=np.float32), (N, N)), -9999.0)

print("  range of bed = [%.2f, %.2f]" % (bed.min(), bed.max()))
print("  range of thk = [%.2f, %.2f]" % (thk.min(), thk.max()))
print("  range of msk = [%.2f, %.2f]" % (msk.min(), msk.max()))

if showgeom:
    print("showing fields ...")

    fig = plt.figure(1)
    ax = plt.imshow(msk)
    fig.colorbar(ax)

    fig = plt.figure(2)
    ax = plt.imshow(bed)
    fig.colorbar(ax)

    fig = plt.figure(3)
    ax = plt.imshow(thk)
    fig.colorbar(ax)

    plt.show()

    sys.exit(0)

#err = abs(bed[msk<0.5] + thk[msk<0.5] - srf[msk<0.5])
# print err.max()

outname = 'ant1kmgeom.nc'

print("writing NetCDF file '%s' ..." % outname)
try:
    nc = PNC(outname, 'w', format='NETCDF3_CLASSIC')
except:
    print("can't open file %s for writing" % outname)
    exit(1)

print("  writing x,y ...")
dx = 1000.0
dy = 1000.0
x = np.linspace(0.0, (N - 1) * dx, N)
y = np.linspace(0.0, (N - 1) * dy, N)
nc.create_dimensions(x, y, time_dependent=False)

print("  writing topg ...")
nc.define_2d_field("topg", time_dependent=False,
                   attrs={"long_name": "elevation of bedrock",
                          "valid_range": (-9000.0, 9000.0),
                          "standard_name": "bedrock_altitude",
                          "units": "meters"})
nc.write_2d_field("topg", bed)

print("  writing thk ...")
nc.define_2d_field("thk", time_dependent=False,
                   attrs={"long_name": "thickness of ice sheet or ice shelf",
                          "valid_range": (0.0, 9000.0),
                          "standard_name": "land_ice_thickness",
                          "units": "meters"})
nc.write_2d_field("thk", thk)

nc.close()
print("done")
