#!/usr/bin/env python

import PISM
from PISM.testing import shallow_grid

ctx = PISM.Context()
config = ctx.config

config.set_flag("stress_balance.calving_front_stress_bc", True)

grid = shallow_grid(Mx=5, My=5, Lx=10e3, Ly=10e3)

geometry = PISM.Geometry(grid)

# grid center:
c = 2
# ice thickness:
H = 1000.0

with PISM.vec.Access(nocomm=[geometry.ice_thickness]):
    for (i, j) in grid.points():
        if abs(i - c) < c and abs(j - c) < c:
            geometry.ice_thickness[i, j] = H
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
geometry.ice_thickness.write(f)
geometry.bed_elevation.write(f)
geometry.cell_type.write(f)
enthalpy.write(f)
tauc.write(f)
f.close()
