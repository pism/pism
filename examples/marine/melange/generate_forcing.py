#!/usr/bin/env python3

from PISMNC import PISMDataset as NC
import numpy as np
import sys

output_filename = sys.argv[1]

with NC(output_filename, 'w') as nc:

    # generate sea level data
    P0 = 1e5
    T = 400.0
    time = np.linspace(0, T, 401)
    delta_MBP = P0 * (1.0 + np.sin(4 * np.pi * time / T))

    dt = time[1] - time[0]
    time_bounds = np.zeros((len(time), 2))
    for k in range(len(time)):
        time_bounds[k, 0] = time[k] - 0.5*dt
        time_bounds[k, 1] = time[k] + 0.5*dt

    nc.create_time(use_bounds=True, length=len(time), units="common_years since 1-1-1")
    nc.write_timeseries("time", time)
    nc.write_timeseries("time_bounds", time_bounds)

    nc.define_timeseries("delta_MBP", attrs={"units": "Pa"})
    nc.write_timeseries("delta_MBP", delta_MBP)
