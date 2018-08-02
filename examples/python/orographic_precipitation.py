#!/usr/bin/env python

import PISM
from math import cos, pi
import numpy as np

ctx = PISM.Context()
config = ctx.config

def gaussian_bump(xmin, xmax, ymin, ymax, dx, dy, h_max=500.0,
                  x0=-25e3, y0=0.0, sigma_x=15e3, sigma_y=15e3):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    # Reproduce Fig 4c in SB2004
    x = np.arange(xmin, xmax, dx)
    y = np.arange(ymin, ymax, dy)
    X, Y = np.meshgrid(x, y)
    Orography = h_max * np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                                 ((Y - y0)**2 / (2 * sigma_y**2))))
    return X, Y, Orography


def initialize_thickness(thickness, H):
    grid = thickness.grid()
    with PISM.vec.Access(nocomm=[thickness]):
        for (i, j) in grid.points():
            thickness[i, j] = orography[i, j]

def allocate(grid):
    H = PISM.model.createIceThicknessVec(grid)
    bed = PISM.model.createBedrockElevationVec(grid)
    sea_level = PISM.IceModelVec2S(grid, "sea_level", PISM.WITHOUT_GHOSTS)

    return H, bed, sea_level

def create_grid():
    P = PISM.GridParameters(config)
    P.horizontal_size_from_options()
    P.horizontal_extent_from_options()
    P.vertical_grid_from_options(config)
    P.ownership_ranges_from_options(ctx.size)

    return PISM.IceGrid(ctx.ctx, P)

def run(plot, pause, save):

    # set grid defaults
    config.set_double("grid.Mx", 400)
    config.set_double("grid.My", 400)

    config.set_double("grid.Lx", 300e3)
    config.set_double("grid.Ly", 200e3)

    config.set_double("grid.Mz", 2)
    config.set_double("grid.Lz", 1000)

    config.set_string("atmosphere.orographic_precipitation.file", "~/pism-olympics/data_sets/climate_forcing/ltop_climate_olympics_1000m_dir_220_kg_m-2_yr-1.nc")
    grid = create_grid()

    thickness, bed, sea_level = allocate(grid)

    # set initial geometry
    bed.set(0.0)
    thickness.set(0.0)
    sea_level.set(0.0)
    initialize_thickness(thickness, orography)
    g = PISM.Geometry(grid)
    op = PISM.AtmosphereOrographicPrecipitation(grid)
    op.init(g)
    op.update(g, 0.1, 0.1)

if __name__ == "__main__":

    x_min   = -100e3
    x_max   = 200e3
    y_min   = -150e3
    y_max   = 150e3
    dx      = 750
    dy      = 750
    x0      = -25e3
    y0      = 0.0
    sigma_x = 15e3
    sigma_y = 15e3

    _, _, orography = gaussian_bump(x_min, x_max, y_min, y_max, dx, dy,
                                   x0=x0, y0=y0, sigma_x=sigma_x, sigma_y=sigma_y)

    plot = PISM.OptionBool("-plot", "Plot bed elevation and uplift.")
    save = PISM.OptionBool("-save", "Save final states of the bed elevation and uplift.")
    pause = PISM.OptionBool("-pause", "Pause for 5 seconds to look at runtime 2D plots.")

    run(plot, pause, save)

