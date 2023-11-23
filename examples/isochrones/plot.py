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

def plot_isochrones(ax, f, index=-1):

    H_max = 4000
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
    ax.plot(x, thk, '.', color='black')

    ax.set_ylim(ymin=-100, ymax=H_max)
    ax.set_xlim(xmin=0, xmax=600)
    ax.set_xlabel("distance from the center, km")
    ax.set_ylabel("elevation, m")
    ax.legend()
    ax.set_title(f"EISMINT-II experiment A\nTime: {int(time)} years")

def animate(fig, ax, input_files, output_file):

    writer = FFMpegWriter(fps=20)

    with writer.saving(fig, output_file, dpi=75):
        for filename in input_files:
            print(f"Processing {filename}...")
            with NC.Dataset(filename, 'r') as f:
                N = len(f.variables['time'])
                for k in range(N):
                    print(f"Frame {k}...")
                    ax.clear()
                    plot_isochrones(ax, f, index=k)
                    writer.grab_frame()

def plot_final_frame(fig, ax, input_files, model_file, output_file):
    H_max = 3500
    with NC.Dataset(input_files[-1], 'r') as f:
        plot_isochrones(ax, f, index=-1)

    with NC.Dataset(model_file, 'r') as f:
        z = f.variables['z'][:]
        ax.hlines(z, xmin=0, xmax=600,
                  linestyles='solid', colors="gray", label="PISM's vertical grid", linewidths=0.5)

    ax.legend()
    ax.set_ylim(ymin=-100, ymax=H_max)
    ax.set_xlim(xmin=0, xmax=600)
    fig.savefig(output_file, dpi=200)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog='plot.py',
        description="Visualizes isochrones using PISM's diagnostic 'isochrone_depth'")
    parser.add_argument('files', metavar='filename', type=str, nargs='+',
                        help='files to process')
    parser.add_argument('-o', dest='output', type=str, nargs=1,
                        help='output file name', required=True)
    parser.add_argument('-p', dest='model_output', type=str, nargs=1,
                        help='model state file')
    parser.add_argument('-f', dest='final', action="store_true",
                        help='plot the final frame from "files"')

    args = parser.parse_args()

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)

    if args.final:
        plot_final_frame(fig, ax, args.files, args.model_output[0], args.output[0])
    else:
        animate(fig, ax, args.files, args.output[0])
