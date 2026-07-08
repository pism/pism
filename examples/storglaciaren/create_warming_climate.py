#!/usr/bin/env python3
# Copyright (C) 2017, 2020, 2024 Andy Aschwanden

import os
import sys
import numpy as np
import time
import xarray as xr
from argparse import ArgumentParser


# Set up the option parser
parser = ArgumentParser()
parser.description = "Create climate forcing for a warming climate"
parser.add_argument("FILE", nargs='*')
parser.add_argument("-T_max", dest="T_max", type=float,
                    help="Maximum temperature", default=1)
parser.add_argument("-t_max", dest="t_max", type=float,
                    help="lower time bound for maximum temperature", default=100)
parser.add_argument("-amplitude", dest="amplitude", type=float,
                    help="Amplitde of seasonal cycle.", default=12)


options = parser.parse_args()
args = options.FILE
start = 0
end = 1000
step = 1./12.
amplitude = options.amplitude
t_max = options.t_max
T_max = options.T_max
bnds_interval_since_refdate = np.linspace(start, end, end * 12 + 1)
time_interval_since_refdate = (bnds_interval_since_refdate[0:-1] +
                               np.diff(bnds_interval_since_refdate) / 2)

infile = args[0]

# Build delta_T time series
T_0 = 0.
temp = np.zeros_like(time_interval_since_refdate) + T_max
temp[0:int(t_max/step)] = np.linspace(T_0, T_max, int(t_max / step))
temp[:] += -np.cos(time_interval_since_refdate * 2 * np.pi) * amplitude

# Build the time-bounds array (Nt, 2)
time_bnds = np.column_stack([bnds_interval_since_refdate[:-1],
                             bnds_interval_since_refdate[1:]])

ds = xr.Dataset(
    coords={
        "time": ("time", time_interval_since_refdate.astype("f8"), {
            "bounds": "time_bnds",
            "units": "years since 1-1-1",
            "calendar": "365_day",
            "standard_name": "time",
            "axis": "T",
        }),
    },
    data_vars={
        "time_bnds": (("time", "nb2"), time_bnds.astype("f8")),
        "delta_T":   (("time",), temp.astype("f4"), {"units": "kelvin"}),
    },
    attrs={
        "history": " ".join([time.ctime(), ":", os.path.basename(__file__),
                              " ".join([str(x) for x in args])]),
        "Conventions": "CF 1.6",
    },
)
ds.to_netcdf(infile, mode="w", unlimited_dims=["time"])
