#!/usr/bin/env python3

import argparse

import numpy as np
from PISMNC import PISMDataset as NC

parser = argparse.ArgumentParser()
parser.description = "Generates an input file for the 'flood' experiment."

parser.add_argument("-M", dest="M", type=int, help="grid size", default=301)
parser.add_argument("-L", dest="L", type=float, help="domain size, meters", default=1e6)
parser.add_argument("-o", dest="output", help="output file name", default="flood.nc")
options = parser.parse_args()
L = options.L
M = options.M

x = np.linspace(-L, L, M)
y = np.linspace(-L, L, M)

xx, yy = np.meshgrid(x, y)

z = np.zeros_like(xx)

thk = np.zeros_like(z)
for phi in np.linspace(0, 2 * np.pi, 4):
    R = 0.5 * L
    r = 0.2 * L
    center_x = R * np.cos(phi)
    center_y = R * np.sin(phi)
    T = phi / (2.0 * np.pi) * 1500.0
    thk += (T / r) * np.sqrt(
        np.maximum(r**2 - (xx - center_x) ** 2 - (yy - center_y) ** 2, 0)
    )

try:
    nc = NC(options.output, "w")
except:
    nc = NC(options.output, "a")

try:
    nc.create_dimensions(x, y, time_dependent=False)

    nc.define_2d_field("topg", attrs={"units": "m", "long_name": "bedrock topography"})
    nc.define_2d_field("thk", attrs={"units": "m", "long_name": "ice thickness"})

    nc.define_2d_field("climatic_mass_balance", attrs={"units": "kg m-2 year-1"})
    nc.define_2d_field("ice_surface_temp", attrs={"units": "Kelvin"})
except:
    pass

nc.write("topg", z)
nc.write("thk", thk)
nc.write("climatic_mass_balance", np.zeros_like(xx))
nc.write("ice_surface_temp", np.zeros_like(xx) + 273.15 - 30.0)

# generate sea level data
time = np.linspace(0, 1000, 1001)
sea_level = np.linspace(0, 700, 1001)

nc.createDimension("time")
nc.define_timeseries("time", attrs={"units": "years since 1-1-1"})
nc.write_timeseries("time", time)

nc.define_timeseries("delta_SL", attrs={"units": "meters"})
nc.write_timeseries("delta_SL", sea_level)

nc.close()
