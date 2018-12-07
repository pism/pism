#!/usr/bin/env python

import netCDF4
import numpy as np

M = 3                           # grid size
L = 1e5                         # domain size
b0 = 0.0                        # bed elevation
H0 = 200.0                      # ice thickness
T_mean_annual = 268.15        # mean annual temperature, Kelvin
T_amplitude = 6             # surface temperature aplitude, Kelvin
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


f = netCDF4.Dataset("input.nc", "w")

f.createDimension("x", Mx)
x = f.createVariable("x", np.float64, ("x",))
x.units = "meters"
x[:] = np.linspace(-Lx, Lx, Mx)

f.createDimension("y", My)
y = f.createVariable("y", np.float64, ("y",))
y.units = "meters"
y[:] = np.linspace(-Ly, Ly, My)

f.createDimension("time", n_records)

time = np.linspace(0, 1, n_records) * seconds_per_year

t = f.createVariable("time", np.float64, ("time",))
t.units = "seconds since 1-1-1"
t[:] = time

H = f.createVariable("thk", np.float64, ("y", "x"))
H.units = "meters"
H.long_name = "ice thickness"
H[:] = H0

b = f.createVariable("topg", np.float64, ("y", "x"))
b.units = "meters"
b.long_name = "bed elevation"
b[:] = b0

T_s = f.createVariable("ice_surface_temp", np.float64, ("y", "x"))
T_s.units = "Kelvin"
T_s[:] = T_mean_annual

M = f.createVariable("climatic_mass_balance", np.float64, ("y", "x"))
M.units = "kg m-2 s-1"
M[:] = M0

dT = f.createVariable("delta_T", np.float64, ("time",))
dT.units = "Kelvin"
dT[:] = T_surface(time, 0.0, T_amplitude, summer_peak_day)

f.close()
