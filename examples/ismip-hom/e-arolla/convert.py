#!/usr/bin/env python3

# Converts tab-delimited inputs into a form that can be used by PISM.

import xarray as xr
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])

x       = data[:, 0]
bed     = data[:, 1]
surface = data[:, 2]
sliding = data[:, 3]

dx = x[1] - x[0]
high_beta = 1e21


def convert(output_filename, no_slip):
    if no_slip:
        tauc = np.full((3, len(x)), high_beta)
    else:
        tauc = np.tile((1.0 - sliding) * high_beta, (3, 1))

    ds = xr.Dataset(
        coords={
            "x": ("x", x.astype("f8"), {"units": "m"}),
            "y": ("y", np.array([-dx, 0.0, dx], dtype="f8"), {"units": "m"}),
        },
        data_vars={
            "topg": (("y", "x"), np.tile(bed, (3, 1)).astype("f8"), {"units": "m"}),
            "thk":  (("y", "x"),
                     np.tile(np.maximum(surface - bed, 0.0), (3, 1)).astype("f8"),
                     {"units": "m"}),
            "tauc": (("y", "x"), tauc.astype("f8"), {"units": "Pa"}),
        },
    )
    ds.to_netcdf(output_filename, mode="w")


convert("arolla100-sliding.nc", no_slip=False)
convert("arolla100-no-slip.nc", no_slip=True)
