#!/usr/bin/env python3

import sys
import time
import xarray as xr
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

    coords = {
        "x": (
            ["x"],
            x,
            {
                "units": "m",
                "axis": "X",
                "standard_name": "projection_x_coordinate",
                "long_name": "x-coordinate in projected coordinate system",
            },
        ),
        "y": (
            ["y"],
            y,
            {
                "units": "m",
                "axis": "Y",
                "standard_name": "projection_y_coordinate",
                "long_name": "y-coordinate in projected coordinate system",
            },
        ),
    }

    ds = xr.Dataset(
        {
            "topg": xr.DataArray(
                data=topg,
                dims=["y", "x"],
                coords=coords,
                attrs={
                    "standard_name": "bedrock_altitude",
                    "units": "m",
                },
            ),
            "thk": xr.DataArray(
                data=thk,
                dims=["y", "x"],
                coords=coords,
                attrs={
                    "standard_name": "land_ice_thickness",
                    "units": "m",
                },
            ),
        },
        attrs={"Conventions": "CF-1.7"},
    )
    ds.to_netcdf(filename)




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

    coords = {
        "time": (["time"], time, {"units": "common_years since 0001-1-1",
                                  "axis": "T",
                                  "calendar": "365_day",
                                  "bounds": "time_bounds",
                                  "_FillValue": False}),
    }


    ds = xr.Dataset(
        {
            "delta_T": xr.DataArray(
                data=delta_T,
                dims=["time"],
                coords=coords,
                attrs={
                    "standard_name": "land_ice_temperature_at_firn_base",
                    "units": "K",
                    "long_name": "Temperature (variation from present)",
                },
            ),
        },
        attrs={"Conventions": "CF-1.7"},
    )
    ds.to_netcdf(filename)



if __name__ == "__main__":
    fill_value = np.nan

    prefix = 'slab'
    ncfile = prefix + '.nc'
    dTfile = prefix + '_dT.nc'

    historystr = time.asctime() + ': ' + ' '.join(sys.argv) + '\n'

    generate_input(ncfile)
    print('Wrote PISM input file %s.' % ncfile)

    generate_forcing(dTfile)
    print('Wrote PISM surface temperature forcing file %s.' % dTfile)
