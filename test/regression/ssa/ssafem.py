#!/usr/bin/env python

import PISM
from PISM.testing import shallow_grid

ctx = PISM.Context()
config = ctx.config

config.set_flag("stress_balance.calving_front_stress_bc", True)

grid = shallow_grid(Mx=7, My=7, Lx=10e3, Ly=10e3)

geometry = PISM.Geometry(grid)

c = 3
H = 1000.0
with PISM.vec.Access(nocomm=[geometry.ice_thickness]):
    for (i, j) in grid.points():
        if i == c and j == c:
            geometry.ice_thickness[i, j] = 1 * H
        elif abs(i - c) < c - 1 and abs(j - c) < c - 1:
            geometry.ice_thickness[i, j] = 1 * H
        elif (i == c and abs(j - c) < c) or (j == c and abs(i - c) < c):
            geometry.ice_thickness[i, j] = 1 * H
        else:
            geometry.ice_thickness[i, j] = 0.0

geometry.bed_elevation.set(0.0)
geometry.sea_level_elevation.set(0.0)
geometry.ice_area_specific_volume.set(0.0)
geometry.ensure_consistency(0.0)

tauc = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
tauc.set(1e6)

enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITH_GHOSTS)
T = 260.0
p = ctx.enthalpy_converter.pressure(1000.0)
enthalpy.set(ctx.enthalpy_converter.enthalpy(T, 0.0, p))

inputs = PISM.StressBalanceInputs()

inputs.geometry = geometry
inputs.basal_melt_rate = None
inputs.melange_back_pressure = None
inputs.basal_yield_stress = tauc
inputs.enthalpy = enthalpy
inputs.age = None

ssa = PISM.SSAFEM(grid)

ssa.init()

ssa.update(inputs, full_update=True)

f = PISM.util.prepare_output(config.get_string("output.file_name"))

ssa.velocity().write(f)
ssa.driving_stress().write(f)
geometry.ice_thickness.write(f)
geometry.bed_elevation.write(f)
geometry.cell_type.write(f)
enthalpy.write(f)
f.close()
