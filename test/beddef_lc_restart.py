#!/usr/bin/env python

""" Sets up and runs the Lingle-Clark bed deformation model.

Compares results of

- Two 1000 year steps
- One 1000 year step, followed by
  * saving the model state and
  * re-initializing
  * one more 1000 year step

Used as a regression test for PISM.LingleClark.
"""

import PISM
from PISM.util import convert
import pylab as plt
import numpy as np
import os

ctx = PISM.Context()

# disc load parameters
disc_radius = convert(1000, "km", "m")
disc_thickness = 1000.0         # meters
# domain size
Lx = 2 * disc_radius
Ly = Lx
N = 101

dt = convert(1000.0, "years", "seconds")


def allocate_storage(grid):
    ice_thickness = PISM.IceModelVec2S(grid, "thk", PISM.WITHOUT_GHOSTS)
    ice_thickness.metadata().set_string("standard_name", "land_ice_thickness")
    bed = PISM.IceModelVec2S(grid, "topg", PISM.WITHOUT_GHOSTS)
    bed_uplift = PISM.IceModelVec2S(grid, "uplift", PISM.WITHOUT_GHOSTS)
    sea_level = PISM.IceModelVec2S(grid, "sea_level", PISM.WITHOUT_GHOSTS)

    return ice_thickness, bed, bed_uplift, sea_level


def add_disc_load(ice_thickness, radius, thickness):
    "Add a disc load with a given radius and thickness."
    grid = ice_thickness.grid()

    with PISM.vec.Access(nocomm=ice_thickness):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            if r <= disc_radius:
                ice_thickness[i, j] = disc_thickness


def run(dt, restart=False):
    "Run the model for 1 time step, stop, save model state, restart, do 1 more step."

    grid = PISM.IceGrid.Shallow(ctx.ctx, Lx, Ly, 0, 0, N, N,
                                PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    model = PISM.LingleClark(grid)

    ice_thickness, bed, bed_uplift, sea_level = allocate_storage(grid)
    grid.variables().add(ice_thickness)

    # start with a flat bed, no ice, and no uplift
    bed.set(0.0)
    bed_uplift.set(0.0)
    ice_thickness.set(0.0)
    sea_level.set(0.0)

    # initialize the model
    model.bootstrap(bed, bed_uplift, ice_thickness, sea_level)

    # add the disc load
    add_disc_load(ice_thickness, disc_radius, disc_thickness)

    # do 1 step
    model.step(ice_thickness, sea_level, dt)

    if restart:
        # save the model state
        filename = "lingle_clark_model_state.nc"
        try:
            PISM.util.prepare_output(filename)
            f = PISM.PIO(grid.com, "netcdf3", filename, PISM.PISM_READWRITE)
            model.write_model_state(f)
            f.close()

            # create a new model
            del model
            model = PISM.LingleClark(grid)

            # initialize
            PISM.PETSc.Options().setValue("-i", filename)
            options = PISM.process_input_options(grid.com, ctx.config)
            model.init(options, ice_thickness, sea_level)
        finally:
            os.remove(filename)

    # do 1 more step
    model.step(ice_thickness, sea_level, dt)

    return model


def compare_vec(v1, v2):
    "Compare two vecs."
    try:
        np.testing.assert_equal(v1.numpy(), v2.numpy())
        print("{} is the same".format(v1.get_name()))
    except:
        plt.figure(figsize=(10, 10))
        diff = v1.numpy() - v2.numpy()
        max_diff = np.max(np.fabs(diff))
        m = plt.imshow(diff, origin="lower")
        plt.colorbar(m)
        plt.title("{}, max. difference {}".format(v1.get_name(), max_diff))


def compare(model1, model2):
    "Compare two models"
    compare_vec(model1.bed_elevation(), model2.bed_elevation())
    compare_vec(model1.uplift(), model2.uplift())
    compare_vec(model1.total_displacement(), model2.total_displacement())
    compare_vec(model1.viscous_displacement(), model2.viscous_displacement())
    compare_vec(model1.relief(), model2.relief())


def lingle_clark_restart_test():
    "Compare straight and re-started runs."
    compare(run(dt),
            run(dt, restart=True))
