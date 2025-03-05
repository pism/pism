#!/usr/bin/env python3

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
import numpy as np
import os

ctx = PISM.Context()

# silence initialization messages
ctx.log.set_threshold(1)

# disc load parameters
disc_radius = convert(1000, "km", "m")
disc_thickness = 1000.0         # meters
# domain size
Lx = 2 * disc_radius
Ly = Lx
N = 101

ctx.config.set_number("bed_deformation.lc.grid_size_factor", 2)

dt = convert(1000.0, "years", "seconds")

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

    grid = PISM.Grid.Shallow(ctx.ctx, Lx, Ly, 0, 0, N, N,
                                PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    model = PISM.LingleClark(grid)

    geometry = PISM.Geometry(grid)

    bed_uplift = PISM.Scalar(grid, "uplift")

    # start with a flat bed, no ice, and no uplift
    geometry.bed_elevation.set(0.0)
    geometry.ice_thickness.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ensure_consistency(0.0)

    bed_uplift.set(0.0)

    # initialize the model
    model.bootstrap(geometry.bed_elevation, bed_uplift, geometry.ice_thickness, geometry.sea_level_elevation)

    # add the disc load
    add_disc_load(geometry.ice_thickness, disc_radius, disc_thickness)

    # do 1 step
    model.step(geometry.ice_thickness, dt)

    if restart:
        # save the model state
        filename = "lingle_clark_model_state.nc"
        try:
            PISM.util.prepare_output(filename)
            f = PISM.File(grid.com, filename, PISM.PISM_NETCDF3, PISM.PISM_READWRITE)
            model.write_model_state(f)
            f.close()

            # create a new model
            del model
            model = PISM.LingleClark(grid)

            # initialize
            ctx.config.set_string("input.file", filename)
            options = PISM.process_input_options(grid.com, ctx.config)
            model.init(options, geometry.ice_thickness, geometry.sea_level_elevation)
        finally:
            os.remove(filename)

    # do 1 more step
    model.step(geometry.ice_thickness, dt)

    return model

def compare_vec(v1, v2):
    "Compare two vecs."
    print("Comparing {}".format(v1.get_name()))
    np.testing.assert_equal(v1.to_numpy(), v2.to_numpy())

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
