#!/usr/bin/env python3

"""Simple testing program for Schoof (2003)-type bed smoothing and
roughness- parameterization schemes. Allows comparison of computed
theta to result from Matlab/Octave code exampletheta.m in
src/base/bedroughplay. Also used in PISM software (regression) test.
"""


import PISM
from math import sin, pi
import numpy

ctx = PISM.Context()
config = ctx.config


def grid():
    "Create the bed smoother grid."
    P = PISM.GridParameters(config)

    P.horizontal_size_and_extent_from_options(config)
    P.vertical_grid_from_options(config)
    P.ownership_ranges_from_options(config, ctx.size)

    return PISM.Grid(ctx.ctx, P)


def allocate_storage(grid):
    "Allocate the bed, the smoothed bed, the surface elevation, and theta."
    topg = PISM.Scalar1(grid, "topg")
    topg.metadata(0).long_name("original topography").units("m").standard_name("bedrock_altitude")

    topg_smoothed = PISM.Scalar(grid, "topg_smoothed")
    topg_smoothed.metadata(0).long_name("smoothed topography").units("m").standard_name("bedrock_altitude")

    usurf = PISM.Scalar1(grid, "usurf")
    usurf.metadata(0).long_name("ice surface elevation").units("m").standard_name("surface_altitude")

    theta = PISM.Scalar1(grid, "theta")
    theta.metadata(0).long_name("coefficient theta in Schoof (2003) bed roughness parameterization")

    return (topg, topg_smoothed, usurf, theta)


def set_topg(topg):
    "Initialize the bed topography."
    grid = topg.grid()

    with PISM.vec.Access(comm=[topg]):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            topg[i, j] = (400.0 * sin(2.0 * pi * x / 600.0e3) +
                          100.0 * sin(2.0 * pi * (x + 1.5 * y) / 40.0e3))


def set_usurf(usurf):
    "Initialize the surface elevation."
    usurf.set(1000.0)


def set_config():
    "Set configuration parameters."

    config.set_string("grid.periodicity", "none")
    config.set_string("grid.registration", "corner")

    config.set_number("grid.Mx", 81)
    config.set_number("grid.My", 81)

    config.set_number("grid.Lx", 1200e3)
    config.set_number("grid.Ly", 1200e3)

    config.set_number("stress_balance.sia.Glen_exponent", 3.0)
    config.set_number("stress_balance.sia.bed_smoother.range", 50.0e3)

    PISM.set_config_from_options(ctx.unit_system, config)


def smooth(topg, topg_smoothed, usurf, theta):
    "Smooth the bed topography."
    grid = topg.grid()

    smoother = PISM.BedSmoother(grid)

    smoother.preprocess_bed(topg)

    smoother.theta(usurf, theta)

    topg_smoothed.copy_from(smoother.smoothed_bed())


def run():
    "Run the bed smoother using synthetic geometry."

    set_config()

    topg, topg_smoothed, usurf, theta = allocate_storage(grid())

    set_usurf(usurf)

    set_topg(topg)

    smooth(topg, topg_smoothed, usurf, theta)

    return topg, topg_smoothed, usurf, theta


def bed_smoother_test():
    "Compare the range of topg, topg_smoothed, and theta to stored values"

    topg, topg_smoothed, usurf, theta = run()

    stored_range = {}
    stored_range["topg"] = [-500.0, 500.0]
    stored_range["topg_smoothed"] = [-372.9924735817933, 372.9924735817933]
    stored_range["theta"] = [0.7147300652935706, 0.9884843647808601]

    for f in [topg, topg_smoothed, theta]:
        numpy.testing.assert_almost_equal([PISM.min(f), PISM.max(f)], stored_range[f.get_name()])


if __name__ == "__main__":
    for field in run():
        field.dump("bed_smoother_%s.nc" % field.get_name())
