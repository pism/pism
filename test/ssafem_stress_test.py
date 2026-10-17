#!/usr/bin/env python3

# This script checks if the SSAFEM solver can handle all possible configurations of ice in
# a 4x4 grid to confirm that the bug in the node type computation that lead to solver
# failures is fixed.
#
# Prior to this fix the SSAFEM solver with CFBC "on" could not handle configurations like
#
# o--o--o--o--o
# |  |  |  |  |
# 1--2--3--4--o
# |  |  |  |  |
# 5--6--o--o--o
# |  |  |  |  |
# o--o--o--o--o
#
# where numbered nodes are "icy" and nodes marked with "o" are ice-free because the node
# "4" does not belong to any "icy" element.
#

import numpy as np

import PISM
import PISM.testing

ctx = PISM.Context()
ctx.config.set_flag("stress_balance.calving_front_stress_bc", True)
ctx.log.set_threshold(1)

grid = PISM.testing.shallow_grid(Mx=6, My=6, Lx=1e4, Ly=1e4)

geometry = PISM.Geometry(grid)
geometry.sea_level_elevation.set(0.0)
geometry.bed_elevation.set(0.0)

tauc = PISM.Scalar(grid, "tauc")
tauc.set(1000 * ctx.config.get_number("basal_yield_stress.constant.value"))

enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())

H = 100
E = ctx.enthalpy_converter.enthalpy(260, 0, ctx.enthalpy_converter.pressure(H))
enthalpy.set(E)

inputs = PISM.StressBalanceInputs()
inputs.geometry = geometry
inputs.enthalpy = enthalpy
inputs.basal_yield_stress = tauc


def M(N):
    """Return the mask (icy=1, ice free=0) for the geometric configuration number `N`."""
    return np.array([int(x) for x in np.binary_repr(N, width=16)], dtype=int).reshape(4, 4)


def test_solver(N):
    """Try running the SSAFEM solver with using geometric configuration number `N`."""
    geometry.ice_thickness.set(0.0)
    geometry.ice_thickness.local_part()[3:7, 3:7] = M(N) * H
    geometry.ensure_consistency(0.0)

    try:
        ssa = PISM.SSAFEM(grid)
        ssa.init()
        ssa.update(inputs, True)
        mag = PISM.Scalar(grid, "vel_mag")
        PISM.compute_magnitude(ssa.velocity(), mag)
        return mag
    except RuntimeError:
        return False


def test_node_type(N, prn=True):
    """Print ice and node types masks for config. `N` if `prn` is True, otherwise return them."""
    mask = M(N)
    geometry.ice_thickness.local_part()[3:7, 3:7] = mask
    node = PISM.Scalar(grid, "node_type")
    PISM.compute_node_types(geometry.ice_thickness, 1e-3, node)
    if prn:
        print("Ice thickness mask")
        print(mask)
        print("Node type mask (interior: -1, boundary: 0, exterior: 1)")
        print(node.to_numpy().astype('b')[1:-1, 1:-1])
    else:
        return mask, node.to_numpy().astype('b')


failed = []
for i in range(2**16):
    if test_solver(i):
        print(f"{i}: OK")
    else:
        print(f"{i}: FAILED")
        failed.append(i)

if len(failed) == 0:
    print("All tests passed")
else:
    print("The following configurations failed:")
    for i in failed:
        print(i)
        test_node_type(i, prn=True)
