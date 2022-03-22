#!/usr/bin/env python3

"""This script shows how to use the fracture density model in isolation, i.e. separately
from the rest of PISM.

Given an input file containing ice thickness, bed topography, x and y components of the
SSA velocity (u_bc, v_bc), and the mask vel_bc_mask (ones where fracture density is
set, zero elsewhere) this code will run the fracture density model for a few time steps
and save its outputs to a file.
"""

import PISM
ctx = PISM.Context()

filename = ctx.config.get_string("input.file")

# create the grid using the input file
grid = PISM.IceGrid.FromFile(ctx.ctx, filename, ["thk"], PISM.CELL_CENTER)

# initialize geometric data
geometry = PISM.Geometry(grid)
geometry.ice_thickness.regrid(filename)
geometry.bed_elevation.regrid(filename)
geometry.sea_level_elevation.set(0.0)
geometry.ensure_consistency(0)

# allocate the fracture density model
flow_law = PISM.FlowLawFactory("stress_balance.ssa.", ctx.config, ctx.enthalpy_converter).create()
fracture = PISM.FractureDensity(grid, flow_law)

# initialize it using zero fracture age and density
fracture.initialize()

# read in the velocity field
velocity = PISM.Vector(grid, "_bc")
velocity.set_attrs("", "x-component of the ice velocity", "m s-1", "m s-1", "", 0)
velocity.set_attrs("", "y-component of the ice velocity", "m s-1", "m s-1", "", 1)
velocity.regrid(filename)

# find the longest time step we can take with this velocity field
data = PISM.max_timestep_cfl_2d(geometry.ice_thickness, geometry.cell_type, velocity)
dt = data.dt_max.value()

# read in the BC mask
vel_bc_mask = PISM.IceModelVec2S(grid, "vel_bc_mask")
vel_bc_mask.regrid(filename)

# Set hardness to a constant corresponding to -30C at 0 Pa
hardness = PISM.IceModelVec2S(grid, "hardness")
T = -30
P = 0.0
E = ctx.enthalpy_converter.enthalpy(T + 273.15, 0.0, P)
hardness.set(flow_law.hardness(E, P))

# take a few steps
N = 100                         # arbitrary
for k in range(N):
    fracture.update(dt, geometry, velocity, hardness, vel_bc_mask)

# save results
output_filename = ctx.config.get_string("output.file")

f = PISM.util.prepare_output(output_filename, append_time=True)
fracture.density().write(f)
fracture.age().write(f)
fracture.growth_rate().write(f)
fracture.flow_enhancement().write(f)
fracture.healing_rate().write(f)
fracture.toughness().write(f)
f.close()
