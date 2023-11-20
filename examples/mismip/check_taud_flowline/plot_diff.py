#!/usr/bin/env python3

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from matplotlib import cm

# ICE FRONT

data_path = ""

experiments = [data_path + "result_1a_A1_800_300_1750e3_601_"]
labels = ["old", "new"]

# FIXME resolution dependent!!

with nc.Dataset(experiments[0] + "default.nc", "r") as ncr:
    mask = np.ma.array(ncr.variables["mask"][0, 1, :])
    x = np.ma.array(ncr.variables["x"][:]) / 1000.0  # convert to km

# find the index of the last grounded cell gl_idx:
MASK_FLOATING = 3
gl_idx = 0
for i in range(len(x) - 1):
    if mask[i] < MASK_FLOATING and mask[i + 1] >= MASK_FLOATING:
        gl_idx = i
        break

# MISMIP
secpera = 3.15569259747e7


def plot_thk(axis, case, i):
    varname = "thk"
    axis.axvline(x=(x[gl_idx]), color="grey", linestyle="--", label=None)
    axis.axhline(y=0, color="grey", linestyle="--", label=None)

    run = experiments[0]

    label = case.split("_")[-1].rstrip(".nc")

    def plot(filename, label):
        with nc.Dataset(filename, "r") as ncr:
            thk = np.ma.array(ncr.variables[varname][0, 1, :])
            thk.mask = mask > MASK_FLOATING

        axis.plot(x, thk, linewidth=2, label=label)
        axis.scatter(x, thk)

    if case != "default.nc":
        plot(run + "default.nc", "default")

    plot(run + case, label)

    axis.set_ylim([0, 100])  # m
    axis.set_xlim([1730, 1750])  # km
    axis.legend(loc="lower right")
    axis.grid(True)

    if i == 0:
        axis.set_title("Change in ice thickness (m)")
    if i < 2:
        axis.set_xticklabels([])


def plot_taud(axis, case, i):
    varname = "taud_x"
    axis.axvline(x=(x[gl_idx]), color="grey", linestyle="--", label=None)
    axis.axhline(y=0, color="grey", linestyle="--", label=None)
    for run in experiments:
        with nc.Dataset(run + "default.nc", "r") as ncr:
            refvariable = np.ma.array(ncr.variables[varname][0, 1, :])
            refvariable.mask = mask > MASK_FLOATING

        with nc.Dataset(run + case, "r") as ncr:
            variable = np.ma.array(ncr.variables[varname][0, 1, :])
            variable.mask = mask > MASK_FLOATING

        axis.plot(x, variable - refvariable, linewidth=2)
        axis.scatter(x, variable - refvariable, label=None)

    axis.set_ylim([-700, 700])  # Pa
    axis.set_xlim([1730, 1750])  # km
    axis.grid(True)

    if i == 0:
        axis.set_title("Change in driving stress (Pa)")
    if i < 2:
        axis.set_xticklabels([])


def plot_velocity(axis, case, i):
    varname = "u_ssa"
    axis.axvline(x=x[gl_idx], color="grey", linestyle="--", label=None)
    axis.axhline(y=0, color="grey", linestyle="--", label=None)
    for run in experiments:
        with nc.Dataset(run + "default.nc", "r") as ncr:
            refvariable = np.ma.array(ncr.variables[varname][0, 1, :])

        with nc.Dataset(run + case, "r") as ncr:
            variable = np.ma.array(ncr.variables[varname][0, 1, :])
            variable.mask = mask > MASK_FLOATING

        axis.plot(x, (variable - refvariable) * secpera, linewidth=2)
        axis.scatter(x, (variable - refvariable) * secpera, label=None)

    axis.set_xlim([1730, 1750])  # km
    axis.grid(True)

    if i == 0:
        axis.set_title("Change in velocity (m/a)")
    if i < 2:
        axis.set_xticklabels([])


fig, axes = plt.subplots(3, 3, sharex=False, sharey=False, figsize=(10, 8))
plt.subplots_adjust(wspace=0.25, hspace=0.1)

for i, case in enumerate(["p1.nc", "p2.nc", "p3.nc"]):
    # plot ice thickness
    plot_thk(axes[i, 0], case, i)

    # plot changes is taud
    plot_taud(axes[i, 1], case, i)

    # plot changes in ice velocity
    plot_velocity(axes[i, 2], case, i)

for k in [0, 1, 2]:
    axes[2, k].set_xlabel("x (km)")

plt.show()
