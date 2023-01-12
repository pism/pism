#!/usr/bin/env python3

"""Creates an animation visualizing isochrones.

Uses PISM's 'isochrone_depth' and 'thk' diagnostics.

Requires `ffmpeg`.

"""

import netCDF4 as NC
import numpy as np
from matplotlib import pyplot as plt
import sys

import matplotlib
matplotlib.use("Agg")

from matplotlib.animation import FFMpegWriter

def plot(ax, f, index=-1):

    Mx = len(f.variables['x'])
    My = len(f.variables['y'])

    deposition_time = f.variables['deposition_time'][:2]
    spacing = (deposition_time[1] - deposition_time[0]) / (365 * 86400)

    time = f.variables['time'][index] / (365 * 86400)
    x = f.variables['x'][Mx // 2:] / 1000.0

    d = f.variables['isochrone_depth'][index, My // 2, Mx // 2:]

    thk = f.variables['thk'][index, My // 2, Mx // 2:]

    n_layers = d.shape[1]

    for k in range(n_layers):
        data = d[:, k]
        if np.max(data) == 0:
            continue

        ax.plot(x, thk - data, '--', color='blue', linewidth=0.75,
                label=f"isochrones (spacing: {int(spacing)} years)" if k == 0 else None)

    ax.plot(x, thk, color='black', label="ice surface")

    ax.grid()
    ax.set_ylim(ymin=-100, ymax=3100)
    ax.set_xlim(xmax=600)
    ax.set_xlabel("x, km")
    ax.set_ylabel("elevation, m")
    ax.legend()
    ax.set_title("Isothermal SIA, SMB from EISMINT-II experiment A\n" +
                 f"Time: {int(time)} years")

if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)

    writer = FFMpegWriter(fps=20)

    with writer.saving(fig, output_filename, dpi=75):
        with NC.Dataset(input_filename, 'r') as f:
            N = len(f.variables['time'])
            for k in range(N):
                print(f"Frame {k}...")
                ax.clear()
                plot(ax, f, index=k)
                writer.grab_frame()
