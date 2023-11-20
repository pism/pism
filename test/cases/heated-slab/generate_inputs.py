#!/usr/bin/env python3

import sys
import time

import netCDF4
import numpy as np


def generate_input(filename):
    "Generate a slab-on-a-flat-bed input file for PISM."
    # Initialize a square 3x3 grid
    Mx = 3
    My = Mx
    Lx = 50e3  # m
    Ly = Lx

    x = np.linspace(-Lx, Lx, Mx)
    y = np.linspace(-Ly, Ly, My)
    dx = x[1] - x[0]
    dy = dx

    topg = np.zeros((My, Mx))
    thk = np.zeros((My, Mx)) + 1000.0

    # Write the data:
    nc = netCDF4.Dataset(filename, "w")

    # Create dimensions x and y
    nc.createDimension("x", size=Mx)
    nc.createDimension("y", size=My)

    x_var = nc.createVariable("x", "f4", dimensions=("x",))
    x_var.units = "m"
    x_var.standard_name = "projection_x_coordinate"
    x_var[:] = x

    y_var = nc.createVariable("y", "f4", dimensions=("y",))
    y_var.units = "m"
    y_var.standard_name = "projection_y_coordinate"
    y_var[:] = y

    def def_var(nc, name, units, fillvalue):
        var = nc.createVariable(name, "f", dimensions=("y", "x"), fill_value=fillvalue)
        var.units = units
        return var

    bed_var = def_var(nc, "topg", "m", fill_value)
    bed_var.standard_name = "bedrock_altitude"
    bed_var[:] = topg

    thk_var = def_var(nc, "thk", "m", fill_value)
    thk_var.standard_name = "land_ice_thickness"
    thk_var[:] = thk

    # set global attributes
    nc.Conventions = "CF-1.4"
    setattr(nc, "history", historystr)

    nc.close()


def generate_forcing(filename):
    "Generate surface temperature forcing."
    ts_yr = 0
    te_yr = 300000
    dt_yr = 10

    Mt = int((te_yr - ts_yr) / dt_yr + 1)
    delta_T = np.zeros(Mt)
    time = np.linspace(ts_yr, te_yr, Mt)
    time_bounds = np.zeros((Mt, 2))
    time_bounds[:, 0] = np.r_[-dt_yr, time[:-1]]
    time_bounds[:, 1] = time

    delta_T[np.where((time > 100000))] = 25.0
    delta_T[np.where((time >= 150000))] = 0.0

    # Write the data:
    nc = netCDF4.Dataset(dTfile, "w")

    # Create dimensions x and y
    nc.createDimension("time", size=Mt)

    def def_var(nc, name, units, fillvalue):
        var = nc.createVariable(name, "f", dimensions=("time"), fill_value=fillvalue)
        var.units = units
        return var

    t_var = def_var(nc, "time", "common_years since 1-1-1", fill_value)
    t_var[:] = time
    t_var.bounds = "time_bounds"
    t_var.calendar = "365_day"

    nc.createDimension("nv", size=2)
    tb_var = nc.createVariable("time_bounds", "f", dimensions=("time", "nv"))
    tb_var[:] = time_bounds

    dT_var = def_var(nc, "delta_T", "K", fill_value)
    dT_var.standard_name = "land_ice_temperature_at_firn_base"
    dT_var.long_name = "Temperature (variation from present)"
    dT_var[:] = delta_T

    # set global attributes
    nc.Conventions = "CF-1.4"
    setattr(nc, "history", historystr)

    nc.close()


if __name__ == "__main__":
    fill_value = np.nan

    prefix = "slab"
    ncfile = prefix + ".nc"
    dTfile = prefix + "_dT.nc"

    historystr = time.asctime() + ": " + " ".join(sys.argv) + "\n"

    generate_input(ncfile)
    print("Wrote PISM input file %s." % ncfile)

    generate_forcing(dTfile)
    print("Wrote PISM surface temperature forcing file %s." % dTfile)
