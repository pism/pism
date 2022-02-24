#!/usr/bin/env python

import PISM
from PISM.testing import shallow_grid

ctx = PISM.Context()
config = ctx.config

config.set_flag("stress_balance.calving_front_stress_bc", True)

M = 21
# grid center:
c = (M - 1) // 2
# ice thickness:
H = 1000.0
# length (in grid cells) of peninsulas sticking out of the blob
L = 1
w = c

grid = shallow_grid(Mx=M, My=M, Lx=10e3, Ly=10e3)

geometry = PISM.Geometry(grid)

with PISM.vec.Access(nocomm=[geometry.ice_thickness]):
    for (i, j) in grid.points():
        if abs(i - c) + abs(j - c) < w:
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

try:
    ssa.update(inputs, full_update=True)
    print("succeeded")
except:
    print("failed")

v_mag = PISM.IceModelVec2S(grid, "vel_mag", PISM.WITHOUT_GHOSTS)

v_mag.set_to_magnitude(ssa.velocity())

node_type = PISM.IceModelVec2S(grid, "node_type")

PISM.compute_node_types(geometry.ice_thickness, 1.0, node_type)

def save_results(filename):

    f = PISM.util.prepare_output(filename)

    node_type.write(f)
    v_mag.write(f)
    ssa.velocity().write(f)
    geometry.ice_thickness.write(f)
    geometry.bed_elevation.write(f)
    geometry.cell_type.write(f)
    enthalpy.write(f)
    tauc.write(f)
    f.close()

save_results(config.get_string("output.file"))

def plot():
    import pylab as plt
    import numpy as np

    x = np.array(grid.x())
    dx = grid.dx()
    y = np.array(grid.y())
    dy = grid.dy()

    v = ssa.velocity().numpy()

    plt.figure(1)
    plt.imshow(v[:,:,0])
    plt.title("u")
    plt.colorbar()

    plt.figure(2)
    plt.imshow(v[:,:,1])
    plt.title("v")
    plt.colorbar()

    plt.show()

# plot()
