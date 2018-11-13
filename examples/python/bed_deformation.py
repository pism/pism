#!/usr/bin/env python

import PISM
from PISM.util import convert
from math import cos, pi

# Simple testing program for Lingle & Clark bed deformation model.
# Runs go for 150,000 years on 63.5km grid with 100a time steps and Z=2 in L&C model.
# SCENARIOS:  run 'python bed_deformation.py -scenario N' where N=1,2,3,4 as follows
#    (1) dump ice disc on initially level, non-uplifting land, use only viscous
#        half-space model:
#              include_elastic = FALSE, do_uplift = FALSE, H0 = 1000.0
#        center depth b(0,0) should eventually equilibriate to near
#        -1000 * (910/3300) = -275.76 m
#    (2) dump ice disc on initially level, non-uplifting land, use both viscous
#        half-space model and elastic model
#              include_elastic = TRUE, do_uplift = FALSE, H0 = 1000.0
#    (3) never loaded, initially level, uplifting land, use only viscous
#        half-space model (because elastic model gives no additional when no load):
#              include_elastic = FALSE, do_uplift = TRUE, H0 = 0.0
#    (4) dump ice disc on initially level, uplifting land, use both viscous
#        half-space model and elastic model:
#              include_elastic = TRUE, do_uplift = TRUE, H0 = 1000.0;

ctx = PISM.Context()
config = ctx.config

R0 = 1000e3


def initialize_uplift(uplift):
    "Initialize the uplift field."
    grid = uplift.grid()
    peak_uplift = convert(10, "mm/year", "m/second")
    with PISM.vec.Access(nocomm=[uplift]):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            if r < 1.5 * R0:
                uplift[i, j] = peak_uplift * (cos(pi * (r / (1.5 * R0))) + 1.0) / 2.0
            else:
                uplift[i, j] = 0.0


def initialize_thickness(thickness, H0):
    grid = thickness.grid()
    with PISM.vec.Access(nocomm=[thickness]):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            if r < R0:
                thickness[i, j] = H0
            else:
                thickness[i, j] = 0.0


def allocate(grid):
    H = PISM.model.createIceThicknessVec(grid)
    bed = PISM.model.createBedrockElevationVec(grid)
    uplift = PISM.IceModelVec2S()
    uplift.create(grid, "uplift", PISM.WITHOUT_GHOSTS)
    uplift.set_attrs("internal", "bed uplift", "m / second", "")

    sea_level = PISM.IceModelVec2S(grid, "sea_level", PISM.WITHOUT_GHOSTS)

    return H, bed, uplift, sea_level


def create_grid():
    P = PISM.GridParameters(config)
    P.horizontal_size_from_options()
    P.horizontal_extent_from_options()
    P.vertical_grid_from_options(config)
    P.ownership_ranges_from_options(ctx.size)

    return PISM.IceGrid(ctx.ctx, P)


def run(scenario, plot, pause, save):

    # set grid defaults
    config.set_double("grid.Mx", 193)
    config.set_double("grid.My", 129)

    config.set_double("grid.Lx", 3000e3)
    config.set_double("grid.Ly", 2000e3)

    config.set_double("grid.Mz", 2)
    config.set_double("grid.Lz", 1000)

    scenarios = {"1": (False, False, 1000.0),
                 "2": (True,  False, 1000.0),
                 "3": (False, True,  0.0),
                 "4": (True,  True,  1000.0)}

    elastic, use_uplift, H0 = scenarios[scenario]

    print("Using scenario %s: elastic model = %s, use uplift = %s, H0 = %f m" % (scenario, elastic, use_uplift, H0))

    config.set_boolean("bed_deformation.lc.elastic_model", elastic)

    grid = create_grid()

    thickness, bed, uplift, sea_level = allocate(grid)

    # set initial geometry and uplift
    bed.set(0.0)
    thickness.set(0.0)
    sea_level.set(0.0)
    if use_uplift:
        initialize_uplift(uplift)

    time = ctx.ctx.time()

    time.init(ctx.ctx.log())

    model = PISM.LingleClark(grid)

    model.bootstrap(bed, uplift, thickness, sea_level)

    # now add the disc load
    initialize_thickness(thickness, H0)

    dt = convert(100, "365 day", "seconds")

    # the time-stepping loop
    while time.current() < time.end():
        # don't go past the end of the run
        dt_current = min(dt, time.end() - time.current())

        model.update(thickness, sea_level, time.current(), dt_current)

        if plot:
            model.bed_elevation().view(400)
            model.uplift().view(400)

        print("t = %s years, dt = %s years" % (time.date(), time.convert_time_interval(dt_current, "years")))
        time.step(dt_current)

    print("Reached t = %s years" % time.date())

    if pause:
        print("Pausing for 5 seconds...")
        PISM.PETSc.Sys.sleep(5)

    if save:
        model.bed_elevation().dump("bed_elevation.nc")
        model.uplift().dump("bed_uplift.nc")


if __name__ == "__main__":
    scenario = PISM.OptionKeyword("-scenario", "choose one of 4 scenarios", "1,2,3,4", "1")
    plot = PISM.OptionBool("-plot", "Plot bed elevation and uplift.")
    save = PISM.OptionBool("-save", "Save final states of the bed elevation and uplift.")
    pause = PISM.OptionBool("-pause", "Pause for 5 seconds to look at runtime 2D plots.")

    run(scenario.value(), plot, pause, save)


def scenario1_test():
    "Test if scenario 1 runs"
    run("1", False, False, False)


def scenario3_test():
    "Test if scenario 3 runs"
    run("3", False, False, False)
