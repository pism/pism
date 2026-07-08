#!/usr/bin/env python3

import xarray as xr
import numpy as np

M = 3                           # grid size
L = 1e5                         # domain size
b0 = 0.0                        # bed elevation
H0 = 200.0                      # ice thickness
T_mean_annual = 268.15        # mean annual temperature, kelvin
T_amplitude = 6             # surface temperature aplitude, kelvin
summer_peak_day = 365/2
seconds_per_year = 365 * 86400
M0 = 0.0                        # mass balance
n_records = 53

# use a square grid
Mx = M
My = M
Lx = L
Ly = L


def T_surface(time, mean, amplitude, summer_peak_day):
    "Surface temperature (cosine yearly cycle)"
    day_length = 86400
    summer_peak = summer_peak_day * day_length
    year_length = float(365 * day_length)

    t = np.mod(time - summer_peak, year_length) / year_length

    return mean + amplitude * np.cos(2 * np.pi * t)


x = np.linspace(-Lx, Lx, Mx)
y = np.linspace(-Ly, Ly, My)
time = np.linspace(0, 1, n_records) * seconds_per_year

ds = xr.Dataset(
    coords={
        "x": ("x", x, {"units": "meters"}),
        "y": ("y", y, {"units": "meters"}),
        "time": ("time", time, {"units": "seconds since 1-1-1"}),
    },
    data_vars={
        "thk":  (("y", "x"), np.full((My, Mx), H0),
                 {"units": "meters", "long_name": "ice thickness"}),
        "topg": (("y", "x"), np.full((My, Mx), b0),
                 {"units": "meters", "long_name": "bed elevation"}),
        "ice_surface_temp": (("y", "x"), np.full((My, Mx), T_mean_annual),
                             {"units": "kelvin"}),
        "climatic_mass_balance": (("y", "x"), np.full((My, Mx), M0),
                                  {"units": "kg m^-2 s^-1"}),
        "delta_T": (("time",),
                    T_surface(time, 0.0, T_amplitude, summer_peak_day),
                    {"units": "kelvin"}),
    },
)
ds.to_netcdf("input.nc", mode="w")
