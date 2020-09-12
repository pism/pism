#!/usr/bin/env python3

"""Reads ISMIP-HOM data, interpolates to the same grid, and saves to files using NumPy's
savez_compressed().

"""

import glob
import re
import os
import numpy as np
from scipy.interpolate import interp1d, interp2d, griddata

def vx_bd(filename, xs):
    "Load x-velocity for flowline experiments (B, D) and interpolate onto a given grid."
    data = np.loadtxt(filename)
    x = data[:, 0]
    v = data[:, 1]
    # Fix some data using the fact that velocities are periodic in x:
    if np.isnan(v[0]):
        v[0] = v[-1]
    return interp1d(x, v, fill_value="extrapolate")(xs)

def vx_ac(filename, xs):
    "Load x-velocity for 3D experiments (A, C) and interpolate onto a given grid."
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    v = data[:, 2]

    # Remove points with missing data
    xx = []
    yy = []
    vv = []
    for X, Y, Z in zip(x, y, v):
        if not np.isnan(Z):
            xx.append(X)
            yy.append(Y)
            vv.append(Z)

    ys = np.zeros_like(xs) + 0.25

    # method="linear" and "cubic" cause segfaults
    result = griddata((xx, yy), vv, (xs, ys), method="nearest")

    return result

def sample(files, exp, length_scale, xs, func):
    "Sample x-velocity along the flow for plotting."
    result = {}

    for filename in files:
        pattern = ".*/([a-z0-9]+){exp}{length}(_surf)?\\.txt".format(exp=exp, length=length_scale)
        m = re.match(pattern, filename)

        if m is not None:
            model = m.group(1)
            result[model] = func(filename, xs)
        else:
            pass

    return result

if __name__ == "__main__":
    ismip_prefix = "./ismip_all/"

    files = glob.glob(ismip_prefix + "**/*.txt", recursive=True)

    # 401 is the highest grid resolution in ISMIP-HOM data
    N_samples=401
    xs = np.linspace(0, 1, N_samples)

    for ex in "abcd":
        print("Experiment: ", ex)

        sampler = vx_bd if ex in "bd" else vx_ac

        for length_scale in ["005", "010", "020", "040", "080", "160"]:
            print("  Length scale: ", length_scale)

            data = sample(files, ex, length_scale, xs, sampler)

            print("Models: ", sorted(list(data.keys())))

            data["x"] = xs

            filename = "ismip-hom-{}-{}".format(ex, length_scale)

            np.savez_compressed(filename, **data)
