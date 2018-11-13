#!/usr/bin/env python

import PISM

ctx = PISM.Context()

def create_dummy_grid():
    "Create a dummy grid"
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)

grid = create_dummy_grid()

ice_thickness = PISM.model.createIceThicknessVec(grid)

u = PISM.IceModelVec3()
u.create(grid, "u", PISM.WITHOUT_GHOSTS)

v = PISM.IceModelVec3()
v.create(grid, "v", PISM.WITHOUT_GHOSTS)

w = PISM.IceModelVec3()
w.create(grid, "w", PISM.WITHOUT_GHOSTS)

ice_thickness.set(4000.0)
u.set(0.0)
v.set(0.0)
w.set(0.0)

model = PISM.AgeModel(grid, None)
input_options = PISM.process_input_options(ctx.com, ctx.config)
model.init(input_options)

inputs = PISM.AgeModelInputs(ice_thickness, u, v, w)

dt = PISM.util.convert(1, "years", "seconds")

model.update(0, dt, inputs)

model.age().dump("age.nc")
