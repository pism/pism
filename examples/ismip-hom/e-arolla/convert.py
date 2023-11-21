#!/usr/bin/env python3

# Converts tab-delimited inputs into a form that can be used by PISM.

import netCDF4
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])

x       = data[:,0]
bed     = data[:, 1]
surface = data[:, 2]
sliding = data[:, 3]

dx = x[1] - x[0]
high_beta = 1e22

def convert(output_filename, no_slip):
    try:
        f = netCDF4.Dataset(output_filename, "w")
        f.createDimension("x", len(x))
        f.createDimension("y", 3)

        X = f.createVariable("x", "f8", ("x",))
        X.units = "m"
        X[:] = x
        Y = f.createVariable("y", "f8", ("y",))
        Y.units = "m"
        Y[:] = [-dx, 0, dx]

        topography = f.createVariable("topg", "f8", ("y", "x"))
        topography.units = "m"
        for k in range(3):
            topography[k, :] = bed

        thickness = f.createVariable("thk", "f8", ("y", "x"))
        thickness.units="m"
        for k in range(3):
            thickness[k, :] = np.maximum(surface - bed, 0.0)

        TAUC = f.createVariable("tauc", "f8", ("y", "x"))
        TAUC.units = "Pa"
        if no_slip:
            TAUC[:, :] = high_beta
        else:
            for k in range(3):
                TAUC[k, :] = (1.0 - sliding) * high_beta

    finally:
        f.close()

convert("arolla100-sliding.nc", no_slip=False)
convert("arolla100-no-slip.nc", no_slip=True)
