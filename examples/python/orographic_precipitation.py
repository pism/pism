#!/usr/bin/env python

"""This script attempts to reproduce Figure 4c from 'A linear theory of orographic
precipitation' by Smith and Barstad, 2004.

Run it with "-i filename.nc" to use the bed topography in filename.nc instead.
"""

import numpy as np
import pylab as plt

import PISM

def gaussian_bump_grid():
    "Allocate the grid for the synthetic geometry test."
    x_min   = -100e3
    x_max   = 200e3
    y_min   = -150e3
    y_max   = 150e3
    dx      = 750
    dy      = 750

    x0 = (x_max + x_min) / 2.0
    y0 = (y_max + y_min) / 2.0

    Lx = (x_max - x_min) / 2.0
    Ly = (y_max - y_min) / 2.0
    Mx = int((x_max - x_min) / dx)
    My = int((y_max - y_min) / dy)

    return PISM.IceGrid_Shallow(PISM.Context().ctx,
                                Lx, Ly,
                                x0, y0,
                                Mx, My,
                                PISM.CELL_CORNER, PISM.NOT_PERIODIC)

def gaussian_bump(x, y, h_max=500.0,
                  x0=0.0, y0=0.0, sigma_x=15e3, sigma_y=15e3):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    X, Y = np.meshgrid(x, y)
    surface = h_max * np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                              ((Y - y0)**2 / (2 * sigma_y**2))))
    return surface

def synthetic_geometry(grid, orography):
    config = PISM.Context().config

    # set wind speed and direction
    config.set_number("atmosphere.orographic_precipitation.wind_speed", 15.0)
    config.set_number("atmosphere.orographic_precipitation.wind_direction", 270)

    config.set_number("atmosphere.orographic_precipitation.conversion_time", 1000.0)
    config.set_number("atmosphere.orographic_precipitation.fallout_time", 1000.0)
    config.set_number("atmosphere.orographic_precipitation.water_vapor_scale_height", 2500.0)
    config.set_number("atmosphere.orographic_precipitation.moist_stability_frequency", 0.005)
    config.set_number("atmosphere.orographic_precipitation.reference_density", 7.4e-3)

    # eliminate the effect of the Coriolis force
    config.set_number("atmosphere.orographic_precipitation.coriolis_latitude", 0.0)

    model    = PISM.AtmosphereOrographicPrecipitation(grid, PISM.AtmosphereUniform(grid))
    geometry = PISM.Geometry(grid)

    with PISM.vec.Access(nocomm=geometry.ice_thickness):
        for i,j in grid.points():
            geometry.ice_thickness[i, j] = orography[j, i]

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)

    # compute surface elevation from ice thickness and bed elevation
    geometry.ensure_consistency(0)

    model.init(geometry)
    model.update(geometry, 0, 1)

    # convert from mm/s to mm/hour
    return model.mean_precipitation().numpy() * 3600

def input_file(filename):

    ctx = PISM.Context()

    grid = PISM.IceGrid.FromFile(ctx.ctx, filename, ["topg"], PISM.CELL_CENTER)

    geometry = PISM.Geometry(grid)

    geometry.bed_elevation.regrid(filename)
    geometry.ice_thickness.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ensure_consistency(0.0)

    model = PISM.AtmosphereOrographicPrecipitation(grid, PISM.AtmosphereUniform(grid))
    model.init(geometry)
    model.update(geometry, 0, 1)

    model.mean_precipitation().dump(config.get_string("output.file"))

if __name__ == "__main__":
    ctx = PISM.Context()
    config = ctx.config

    if config.get_string("input.file") != "":
        input_file(config.get_string("input.file"))
    else:
        grid = gaussian_bump_grid()
        x = grid.x()
        y = grid.y()
        orography = gaussian_bump(x, y)

        P = synthetic_geometry(grid, orography)

        levels = np.linspace(0.025, 2.025, 6)

        plt.figure()
        plt.contour(x, y, orography, colors="black", linestyles="dashed")
        cs = plt.contour(x, y, P, levels=levels, linestyles="solid")
        plt.clabel(cs)
        plt.grid()
        plt.xlabel("Distance (m)")
        plt.ylabel("Distance (m)")
        plt.title("Figure 4c from Smith and Barstad,\nA linear theory of orographic precipitation, 2004")
        plt.show()
