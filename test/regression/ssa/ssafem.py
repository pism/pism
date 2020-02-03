#!/usr/bin/env python

import PISM
from PISM.testing import shallow_grid

ctx = PISM.Context()
config = ctx.config

config.set_flag("stress_balance.calving_front_stress_bc", True)

grid = shallow_grid(Mx=3, My=3, Lx=10e3, Ly=10e3)

geometry = PISM.Geometry(grid)

# grid center:
c = 1
# ice thickness:
H = 1000.0
# length (in grid cells) of peninsulas sticking out of the blob
L = 1
w = 2

with PISM.vec.Access(nocomm=[geometry.ice_thickness]):
    for (i, j) in grid.points():
        if i == c and j == c:
            geometry.ice_thickness[i, j] = 1 * H
        elif abs(i - c) < w - L and abs(j - c) < w - L:
            geometry.ice_thickness[i, j] = 1 * H
        elif (i == c and abs(j - c) < w) or (j == c and abs(i - c) < w):
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

try:
    ssa.update(inputs, full_update=True)
    print("succeeded")
except:
    print("failed")

v_mag = PISM.IceModelVec2S(grid, "vel_mag", PISM.WITHOUT_GHOSTS)

v_mag.set_to_magnitude(ssa.velocity())

node_type = PISM.IceModelVec2Int(grid, "node_type", PISM.WITHOUT_GHOSTS)

PISM.compute_node_types(geometry.ice_thickness, 1.0, node_type)

f = PISM.util.prepare_output(config.get_string("output.file_name"))

node_type.write(f)
v_mag.write(f)
ssa.velocity().write(f)
geometry.ice_thickness.write(f)
geometry.bed_elevation.write(f)
geometry.cell_type.write(f)
enthalpy.write(f)
tauc.write(f)
f.close()
