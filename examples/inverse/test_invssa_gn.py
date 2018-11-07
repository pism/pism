#! /usr/bin/env python
#
# Copyright (C) 2012, 2014, 2015, 2016, 2017, 2018 David Maxwell
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import sys
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import os
import math

import PISM


def adjustTauc(mask, tauc):
    """Where ice is floating or land is ice-free, tauc should be adjusted to have some preset default values."""

    grid = mask.grid()
    high_tauc = grid.ctx().config().get_double("basal_yield_stress.ice_free_bedrock")

    with PISM.vec.Access(comm=tauc, nocomm=mask):
        for (i, j) in grid.points():
            if mask.ocean(i, j):
                tauc[i, j] = 0
            elif mask.ice_free(i, j):
                tauc[i, j] = high_tauc


# Main code starts here
if __name__ == "__main__":
    context = PISM.Context()
    config = context.config
    com = context.com
    PISM.set_abort_on_sigint(True)

    append_mode = False

    input_filename = config.get_string("input.file")
    inv_data_filename = PISM.OptionString("-inv_data", "inverse data file", input_filename).value()
    use_tauc_prior = PISM.OptionBool("-inv_use_tauc_prior",
                                     "Use tauc_prior from inverse data file as initial guess.")

    ssarun = PISM.invert.ssa.SSAForwardRunFromInputFile(input_filename, inv_data_filename, 'tauc')
    ssarun.setup()

    vecs = ssarun.modeldata.vecs
    grid = ssarun.grid

    # Determine the prior guess for tauc. This can be one of
    # a) tauc from the input file (default)
    # b) tauc_prior from the inv_datafile if -use_tauc_prior is set
    tauc_prior = PISM.model.createYieldStressVec(grid, 'tauc_prior')
    tauc_prior.set_attrs(
        "diagnostic", "initial guess for (pseudo-plastic) basal yield stress in an inversion", "Pa", "")
    tauc = PISM.model.createYieldStressVec(grid)
    if use_tauc_prior:
        tauc_prior.regrid(inv_data_filename, critical=True)
    else:
        if not PISM.util.fileHasVariable(input_filename, "tauc"):
            PISM.verbPrintf(
                1, com, "Initial guess for tauc is not available as 'tauc' in %s.\nYou can provide an initial guess as 'tauc_prior' using the command line option -use_tauc_prior." % input_filename)
            exit(1)
        tauc.regrid(input_filename, True)
        tauc_prior.copy_from(tauc)

    adjustTauc(vecs.ice_mask, tauc_prior)

    # Convert tauc_prior -> zeta_prior
    zeta = PISM.IceModelVec2S()
    WIDE_STENCIL = int(grid.ctx().config().get_double("grid.max_stencil_width"))
    zeta.create(grid, "", PISM.WITH_GHOSTS, WIDE_STENCIL)
    ssarun.tauc_param.convertFromDesignVariable(tauc_prior, zeta)
    ssarun.ssa.linearize_at(zeta)

    vel_ssa_observed = None
    vel_ssa_observed = PISM.model.create2dVelocityVec(grid, '_ssa_observed', stencil_width=2)
    if PISM.util.fileHasVariable(inv_data_filename, "u_ssa_observed"):
        vel_ssa_observed.regrid(inv_data_filename, True)
    else:
        if not PISM.util.fileHasVariable(inv_data_filename, "u_surface_observed"):
            PISM.verbPrintf(
                1, context.com, "Neither u/v_ssa_observed nor u/v_surface_observed is available in %s.\nAt least one must be specified.\n" % inv_data_filename)
            exit(1)
        vel_surface_observed = PISM.model.create2dVelocityVec(grid, '_surface_observed', stencil_width=2)
        vel_surface_observed.regrid(inv_data_filename, True)

        sia_solver = PISM.SIAFD
        if is_regional:
            sia_solver = PISM.SIAFD_Regional
        vel_sia_observed = PISM.sia.computeSIASurfaceVelocities(modeldata, sia_solver)

        vel_sia_observed.metadata(0).set_name('u_sia_observed')
        vel_sia_observed.metadata(0).set_string('long_name', "x-component of the 'observed' SIA velocities")

        vel_sia_observed.metadata(1).set_name('v_sia_observed')
        vel_sia_observed.metadata(1).set_string('long_name', "y-component of the 'observed' SIA velocities")

        vel_ssa_observed.copy_from(vel_surface_observed)
        vel_ssa_observed.add(-1, vel_sia_observed)

    (designFunctional, stateFunctional) = PISM.invert.ssa.createTikhonovFunctionals(ssarun)
    eta = config.get_double("inverse.tikhonov.penalty_weight")

    solver_gn = PISM.InvSSATikhonovGN(ssarun.ssa, zeta, vel_ssa_observed, eta, designFunctional, stateFunctional)

    seed = PISM.OptionInteger("-inv_seed", "random generator seed")
    if seed.is_set():
        np.random.seed(seed.value() + PISM.Context().rank)

    d1 = PISM.vec.randVectorS(grid, 1)
    d2 = PISM.vec.randVectorS(grid, 1)

    GNd1 = PISM.IceModelVec2S()
    GNd1.create(grid, "", PISM.WITHOUT_GHOSTS)
    GNd2 = PISM.IceModelVec2S()
    GNd2.create(grid, "", PISM.WITHOUT_GHOSTS)

    solver_gn.apply_GN(d1, GNd1)
    solver_gn.apply_GN(d2, GNd2)

    ip1 = d1.get_vec().dot(GNd2.get_vec())
    ip2 = d2.get_vec().dot(GNd1.get_vec())
    PISM.verbPrintf(1, grid.com, "Test of Gauss-Newton symmetry (x^t GN y) vs (y^t GN x)\n")
    PISM.verbPrintf(1, grid.com, "ip1 %.10g ip2 %.10g\n" % (ip1, ip2))
    PISM.verbPrintf(1, grid.com, "relative error %.10g\n" % abs((ip1 - ip2) / ip1))
