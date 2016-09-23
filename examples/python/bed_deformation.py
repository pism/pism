#!/usr/bin/env python

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

scenarios = {"1" : (False, False, 1000.0),
             "2" : (True,  False, 1000.0),
             "3" : (False, True,  0.0),
             "4" : (True,  True,  1000.0)}

R0 = 1000e3
t_final_years = 15e3

dt_years = 100.0
dt = PISM.convert(ctx.unit_system, dt_years, "years", "seconds")

one_year = PISM.convert(ctx.unit_system, 1.0, "years", "seconds")

import PISM
from math import cos, pi

ctx = PISM.Context()

config = ctx.config

# set grid defaults
config.set_double("grid.Mx", 193)
config.set_double("grid.My", 129)

config.set_double("grid.Lx", 3000e3)
config.set_double("grid.Ly", 2000e3)

config.set_double("grid.Mz", 2)
config.set_double("grid.Lz", 1000)

scenario = PISM.OptionKeyword("-scenario", "choose one of 4 scenarios", "1,2,3,4", "1")

elastic, use_uplift, H0 = scenarios[scenario.value()]

config.set_boolean("bed_deformation.lc_elastic_model", elastic)

def initialize_uplift(uplift):
    "Initialize the uplift field."
    peak_uplift = PISM.convert(ctx.unit_system, 10, "mm/year", "m/second")
    with PISM.vec.Access(nocomm=[uplift]):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            if r < 1.5 * R0:
                uplift[i, j] = peak_uplift * (cos(pi * (r / (1.5 * R0))) + 1.0) / 2.0;
            else:
                uplift[i, j] = 0.0

def initialize_thickness(thickness):
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

    return H, bed, uplift

def create_grid():
    P = PISM.GridParameters(config)
    P.horizontal_size_from_options()
    P.horizontal_extent_from_options()
    P.vertical_grid_from_options(config)
    P.ownership_ranges_from_options(ctx.size)

    return PISM.IceGrid(ctx.ctx, P)

grid = create_grid()

thickness, bed, uplift = allocate(grid)

# FIXME: this should not be necessary.
grid.variables().add(thickness)
grid.variables().add(bed)
grid.variables().add(uplift)

# set initial geometry and uplift
bed.set(0.0)
thickness.set(0.0)
if use_uplift:
    initialize_uplift(uplift)

model = PISM.PBLingleClark(grid)
model.init()

# FIXME! This does not work for some reason.
# model.init(bed, uplift, thickness)

# now add the disc load
initialize_thickness(thickness)

# the time stepping loop
n_steps = int(t_final_years / dt_years)

for k in range(n_steps + 1):
    time_years = k * dt_years

    time = one_year * time_years

    model.update(thickness, time, dt)

    print "Time: %f years" % time_years

    model.bed_elevation().view(400)
    model.uplift().view(400)

PISM.PETSc.Sys.sleep(5)

model.bed_elevation().dump("bed_elevation.nc")
model.uplift().dump("bed_uplift.nc")
