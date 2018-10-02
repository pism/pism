#!/usr/bin/env python

import sys
import time
import netCDF4
import numpy as np


def generate_input(filename):
    "Generate a slab-on-a-flat-bed input file for PISM."
    # Initialize a square 3x3 grid
    Mx = 3
    My = Mx
    Lx = 50e3   # m
    Ly = Lx

    x = np.linspace(-Lx, Lx, Mx)
    y = np.linspace(-Ly, Ly, My)
    dx = x[1]-x[0]
    dy = dx

    topg = np.zeros((My, Mx))
    thk = np.zeros((My, Mx)) + 1000.0
    artm = np.zeros((My, Mx)) + 273.15 - 30  # -30 degrees Celsius
    hflx = np.zeros((My, Mx)) + 0.042       # 42 mW/m2
    smb = np.zeros((My, Mx))                # zero

    # Write the data:
    nc = netCDF4.Dataset(filename, "w")

    # Create dimensions x and y
    nc.createDimension("x", size=Mx)
    nc.createDimension("y", size=My)

    x_var = nc.createVariable("x", 'f4', dimensions=("x",))
    x_var.units = "m"
    x_var.standard_name = "projection_x_coordinate"
    x_var[:] = x

    y_var = nc.createVariable("y", 'f4', dimensions=("y",))
    y_var.units = "m"
    y_var.standard_name = "projection_y_coordinate"
    y_var[:] = y

    def def_var(nc, name, units, fillvalue):
        var = nc.createVariable(name, 'f', dimensions=("y", "x"), fill_value=fillvalue)
        var.units = units
        return var

    bed_var = def_var(nc, "topg", "m", fill_value)
    bed_var.standard_name = "bedrock_altitude"
    bed_var[:] = topg

    thk_var = def_var(nc, "thk", "m", fill_value)
    thk_var.standard_name = "land_ice_thickness"
    thk_var[:] = thk

    hflx_var = def_var(nc, "bheatflux", "W m-2", fill_value)
    hflx_var.standard_name = "land_ice_basal_heat_flux"
    hflx_var[:] = hflx

    smb_var = def_var(nc, "climatic_mass_balance", "kg m-2 s-1", fill_value)
    smb_var.standard_name = "land_ice_surface_specific_mass_balance"
    smb_var[:] = smb

    artm_var = def_var(nc, "ice_surface_temp", "K", fill_value)
    artm_var[:] = artm

    # set global attributes
    nc.Conventions = 'CF-1.4'
    setattr(nc, 'history', historystr)

    nc.close()


def generate_forcing(filename):
    "Generate surface temperature forcing."
    ts_yr = 0
    te_yr = 300000
    dt_yr = 10

    Mt = (te_yr - ts_yr) / dt_yr + 1
    delta_T = np.zeros(Mt)
    time = np.linspace(ts_yr, te_yr, Mt)

    delta_T[np.where((time > 100000))] = 25.0
    delta_T[np.where((time >= 150000))] = 0.0

    # Write the data:
    nc = netCDF4.Dataset(dTfile, "w")

    # Create dimensions x and y
    nc.createDimension("time", size=Mt)

    def def_var(nc, name, units, fillvalue):
        var = nc.createVariable(name, 'f', dimensions=("time"), fill_value=fillvalue)
        var.units = units
        return var

    t_var = def_var(nc, "time", "years since since 1-1-1", fill_value)
    t_var[:] = time

    dT_var = def_var(nc, "delta_T", "K", fill_value)
    dT_var.standard_name = "land_ice_temperature_at_firn_base"
    dT_var.long_name = "Temperature (variation from present)"
    dT_var[:] = delta_T

    # set global attributes
    nc.Conventions = 'CF-1.4'
    setattr(nc, 'history', historystr)

    nc.close()


if __name__ == "__main__":
    fill_value = np.nan

    prefix = 'box'
    ncfile = prefix + '.nc'
    dTfile = prefix + '_dT.nc'

    historysep = ' '
    historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'

    generate_input(ncfile)
    print('Wrote PISM input file %s.' % ncfile)

    generate_forcing(dTfile)
    print('Wrote PISM surface temperature forcing file %s.' % dTfile)
