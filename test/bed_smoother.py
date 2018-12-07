#!/usr/bin/env python

"""Simple testing program for Schoof (2003)-type bed smoothing and
roughness- parameterization schemes. Allows comparison of computed
theta to result from Matlab/Octave code exampletheta.m in
src/base/bedroughplay. Also used in PISM software (regression) test.
"""


import PISM
from math import sin, pi

ctx = PISM.Context()
config = ctx.config


def grid():
    "Create the bed smoother grid."
    P = PISM.GridParameters(config)

    P.horizontal_size_from_options()
    P.horizontal_extent_from_options()
    P.vertical_grid_from_options(config)
    P.horizontal_extent_from_options()
    P.ownership_ranges_from_options(ctx.size)

    return PISM.IceGrid(ctx.ctx, P)


def allocate_storage(grid):
    "Allocate the bed, the smoothed bed, the surface elevation, and theta."
    topg = PISM.IceModelVec2S()
    topg.create(grid, "topg", PISM.WITH_GHOSTS, 1)
    topg.set_attrs("internal", "original topography",
                   "m", "bedrock_altitude")

    topg_smoothed = PISM.IceModelVec2S()
    topg_smoothed.create(grid, "topg_smoothed", PISM.WITHOUT_GHOSTS, 1)
    topg_smoothed.set_attrs("internal", "smoothed topography",
                            "m", "bedrock_altitude")

    usurf = PISM.IceModelVec2S()
    usurf.create(grid, "usurf", PISM.WITH_GHOSTS, 1)
    usurf.set_attrs("internal", "ice surface elevation",
                    "m", "surface_altitude")

    theta = PISM.IceModelVec2S()
    theta.create(grid, "theta", PISM.WITH_GHOSTS, 1)
    theta.set_attrs("internal",
                    "coefficient theta in Schoof (2003) bed roughness parameterization",
                    "", "")

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

    config.set_double("grid.Mx", 81)
    config.set_double("grid.My", 81)

    config.set_double("grid.Lx", 1200e3)
    config.set_double("grid.Ly", 1200e3)

    config.set_double("stress_balance.sia.Glen_exponent", 3.0)
    config.set_double("stress_balance.sia.bed_smoother.range", 50.0e3)


def smooth(topg, topg_smoothed, usurf, theta):
    "Smooth the bed topography."
    grid = topg.grid()

    smoother = PISM.BedSmoother(grid, 1)

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

    computed_range = {}
    for f in [topg, topg_smoothed, theta]:
        R = f.range()
        computed_range[f.get_name()] = [R.min, R.max]

    for name in list(stored_range.keys()):
        computed = computed_range[name]
        stored = stored_range[name]

        for k in range(2):
            assert abs(computed[k] - stored[k]) < 1e-16


if __name__ == "__main__":
    for field in run():
        field.dump("bed_smoother_%s.nc" % field.get_name())
