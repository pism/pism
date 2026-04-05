#! /usr/bin/env python3
#
# Copyright (C) 2012, 2014, 2015, 2016, 2017, 2018, 2019 David Maxwell
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

# Compares the computations of the key methods of a forward problem
# with alternative computations (finite differences for derivatives, etc.)
# Takes all the arguments of pismi.py, plus -inv_test with the name of
# one of the tests in the file.  The solver is set with its initial value
# of tauc equal to what it would be for pismi, and then the tests occur.
#
# Tests:
# lin -- the linearized forward map, compared to finite difference approximations.
# lin_transpose -- the transpose of this map, T^*, compared to the formula <Td,r> = <d,T^*r>
# j_design -- the derivative of the SSA residual with respect to the design variable, compared to finite difference approx
# J_design_transpose -- the transpose of this map, J^*, compared to the formula <Jd,r> = <d,J^*r>

import sys
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import os
import math
import PISM.invert.ssa
from PISM.logging import logMessage

import PISM


def view(vec, viewer):
    if (isinstance(vec, PISM.Scalar)):
        v_global = PISM.Scalar(vec.grid(), "")
    else:
        v_global = PISM.Vector(vec.grid(), "")
    v_global.copy_from(vec)
    v_global.get_vec().view(viewer)


def adjustTauc(mask, tauc):
    """Where ice is floating or land is ice-free, tauc should be adjusted to have some preset default values."""

    grid = mask.grid()
    high_tauc = grid.ctx().config().get_number("basal_yield_stress.ice_free_bedrock")

    with PISM.vec.Access(comm=tauc, nocomm=mask):
        for (i, j) in grid.points():
            if mask.ocean(i, j):
                tauc[i, j] = 0
            elif mask.ice_free(i, j):
                tauc[i, j] = high_tauc


def createDesignVec(grid, design_var, name=None, **kwargs):
    if name is None:
        name = design_var
    if design_var == "tauc":
        design_vec = PISM.model.createYieldStressVec(grid, name=name, **kwargs)
    elif design_var == "hardav":
        design_vec = PISM.model.createAveragedHardnessVec(grid, name=name, **kwargs)
    else:
        raise ValueError("Unknown design variable %s" % design_var)
    return design_vec


WIDE_STENCIL = 2


def test_lin(ssarun):
    grid = ssarun.grid

    PISM.verbPrintf(1, grid.com, "\nTest Linearization (Comparison with finite differences):\n")

    S = 250
    d_viewer = PETSc.Viewer().createDraw(title="d", size=S)
    Td_viewer = PETSc.Viewer().createDraw(title="Td", size=S)
    Td_fd_viewer = PETSc.Viewer().createDraw(title="Td_fd", size=S)
    d_Td_viewer = PETSc.Viewer().createDraw(title="d_Td", size=S)

    for (i, j) in grid.points():
        d = PISM.Scalar1(grid, "")
        d.set(0)
        with PISM.vec.Access(comm=d):
            d[i, j] = 1

        ssarun.solver.linearize_at(zeta1)
        u1 = PISM.Vector1(grid, "")
        u1.copy_from(ssarun.solver.solution())

        Td = PISM.Vector1(grid, "")
        ssarun.solver.apply_linearization(d, Td)

        eps = 1e-8
        zeta2 = PISM.Scalar1(grid, "")
        zeta2.copy_from(d)
        zeta2.scale(eps)
        zeta2.add(1, zeta1)
        ssarun.solver.linearize_at(zeta2)
        u2 = PISM.Vector1(grid, "")
        u2.copy_from(ssarun.solver.solution())

        Td_fd = PISM.Vector1(grid, "")
        Td_fd.copy_from(u2)
        Td_fd.add(-1, u1)
        Td_fd.scale(1. / eps)

        d_Td = PISM.Vector1(grid, "")
        d_Td.copy_from(Td_fd)
        d_Td.add(-1, Td)

        n_Td_fd = Td_fd.norm(PETSc.NormType.NORM_2)
        n_Td_l2 = Td.norm(PETSc.NormType.NORM_2)
        n_Td_l1 = Td.norm(PETSc.NormType.NORM_1)
        n_Td_linf = Td.norm(PETSc.NormType.NORM_INFINITY)

        n_d_Td_l2 = d_Td.norm(PETSc.NormType.NORM_2)
        n_d_Td_l1 = d_Td.norm(PETSc.NormType.NORM_1)
        n_d_Td_linf = d_Td.norm(PETSc.NormType.NORM_INFINITY)

        PISM.verbPrintf(1, grid.com, "(i,j)=(%d,%d)\n" % (i, j))
        PISM.verbPrintf(1, grid.com, "apply_linearization(d): l2 norm %.10g; finite difference %.10g\n" %
                        (n_Td_l2, n_Td_fd))

        r_d_l2 = 0
        if n_Td_l2 != 0:
            r_d_l2 = n_d_Td_l2 / n_Td_l2
        r_d_l1 = 0
        if n_Td_l1 != 0:
            r_d_l1 = n_d_Td_l1 / n_Td_l1

        r_d_linf = 0
        if n_Td_linf != 0:
            r_d_linf = n_d_Td_linf / n_Td_linf
        PISM.verbPrintf(1, grid.com, "relative difference: l2 norm %.10g l1 norm %.10g linf norm %.10g\n" %
                        (r_d_l2, r_d_l1, r_d_linf))

        PISM.verbPrintf(1, grid.com, "\n")

        d_global = PISM.Scalar(grid, "")
        d_global.copy_from(d)
        d_global.get_vec().view(d_viewer)

        Td_global = PISM.Vector(grid, "")
        Td_global.copy_from(Td)
        Td_global.get_vec().view(Td_viewer)

        Td_fd_global = PISM.Vector(grid, "")
        Td_fd_global.copy_from(Td_fd)
        Td_fd_global.get_vec().view(Td_fd_viewer)

        d_Td_global = PISM.Vector(grid, "")
        d_Td_global.copy_from(d_Td)
        d_Td_global.get_vec().view(d_Td_viewer)

        PISM.logging.pause()

# ######################################################################################################################
# Jacobian design


def test_j_design(ssarun):
    grid = ssarun.grid

    S = 250
    d_viewer = PETSc.Viewer().createDraw(title="d", size=S)
    drhs_viewer = PETSc.Viewer().createDraw(title="drhs", size=S)
    drhs_fd_viewer = PETSc.Viewer().createDraw(title="drhs_fd", size=S)
    d_drhs_viewer = PETSc.Viewer().createDraw(title="d_drhs", size=S)

    ssarun.solver.linearize_at(zeta1)
    u1 = PISM.Vector1(grid, "")
    u1.copy_from(ssarun.solver.solution())

    for (i, j) in grid.points():
        d = PISM.Scalar1(grid, "")
        d.set(0)
        with PISM.vec.Access(comm=d):
            d[i, j] = 1

        ssarun.solver.linearize_at(zeta1)

        rhs1 = PISM.Vector(grid, "")
        ssarun.solver.assemble_residual(u1, rhs1)

        eps = 1e-8
        zeta2 = PISM.Scalar1(grid, "zeta_prior")
        zeta2.copy_from(d)
        zeta2.scale(eps)
        zeta2.add(1, zeta1)
        ssarun.solver.set_design(zeta2)

        rhs2 = PISM.Vector(grid, "")
        ssarun.solver.assemble_residual(u1, rhs2)

        drhs_fd = PISM.Vector(grid, "")
        drhs_fd.copy_from(rhs2)
        drhs_fd.add(-1, rhs1)
        drhs_fd.scale(1. / eps)

        drhs = PISM.Vector(grid, "")
        ssarun.solver.apply_jacobian_design(u1, d, drhs)

        d_drhs = PISM.Vector(grid, "")

        d_drhs.copy_from(drhs)
        d_drhs.add(-1, drhs_fd)

        n_drhs_fd = drhs_fd.norm(PETSc.NormType.NORM_2)
        n_drhs_l2 = drhs.norm(PETSc.NormType.NORM_2)
        n_drhs_l1 = drhs.norm(PETSc.NormType.NORM_1)
        n_drhs_linf = drhs.norm(PETSc.NormType.NORM_INFINITY)

        n_d_drhs_l2 = d_drhs.norm(PETSc.NormType.NORM_2)
        n_d_drhs_l1 = d_drhs.norm(PETSc.NormType.NORM_1)
        n_d_drhs_linf = d_drhs.norm(PETSc.NormType.NORM_INFINITY)

        PISM.verbPrintf(1, grid.com, "\nTest Jacobian Design (Comparison with finite differences):\n")
        PISM.verbPrintf(1, grid.com, "jacobian_design(d): l2 norm %.10g; finite difference %.10g\n" %
                        (n_drhs_l2, n_drhs_fd))
        if n_drhs_linf == 0:
            PISM.verbPrintf(1, grid.com, "difference: l2 norm %.10g l1 norm %.10g linf norm %.10g\n" %
                            (n_d_drhs_l2, n_d_drhs_l1, n_d_drhs_linf))
        else:
            PISM.verbPrintf(1, grid.com, "relative difference: l2 norm %.10g l1 norm %.10g linf norm %.10g\n" %
                            (n_d_drhs_l2 / n_drhs_l2, n_d_drhs_l1 / n_drhs_l1, n_d_drhs_linf / n_drhs_linf))

        view(d, d_viewer)
        view(drhs, drhs_viewer)
        view(drhs_fd, drhs_fd_viewer)
        view(d_drhs, d_drhs_viewer)

        PISM.logging.pause()


def test_j_design_transpose(ssarun):
    grid = ssarun.grid

    S = 250
    r_viewer = PETSc.Viewer().createDraw(title="r", size=S)
    JStarR_viewer = PETSc.Viewer().createDraw(title="JStarR", size=S)
    JStarR_indirect_viewer = PETSc.Viewer().createDraw(title="JStarR (ind)", size=S)
    d_JStarR_viewer = PETSc.Viewer().createDraw(title="d_JStarR_fd", size=S)

    ssarun.solver.linearize_at(zeta1)
    u = PISM.Vector1(grid, "")
    u.copy_from(ssarun.solver.solution())

    Jd = PISM.Vector(grid, "")

    JStarR = PISM.Scalar(grid, "")

    JStarR_indirect = PISM.Scalar(grid, "")

    for (i, j) in grid.points():

        for k in range(2):

            r = PISM.Vector1(grid, "")
            r.set(0)
            with PISM.vec.Access(comm=r):
                if k == 0:
                    r[i, j].u = 1
                else:
                    r[i, j].v = 1

            ssarun.solver.apply_jacobian_design_transpose(u, r, JStarR)

            r_global = PISM.Vector(grid, "")
            r_global.copy_from(r)

            for (k, l) in grid.points():
                with PISM.vec.Access(nocomm=JStarR_indirect):
                    d = PISM.Scalar1(grid, "")
                    d.set(0)
                    with PISM.vec.Access(comm=d):
                        d[k, l] = 1

                    ssarun.solver.apply_jacobian_design(u, d, Jd)

                    JStarR_indirect[k, l] = Jd.get_vec().dot(r_global.get_vec())

            d_JStarR = PISM.Scalar(grid, "")

            d_JStarR.copy_from(JStarR)
            d_JStarR.add(-1, JStarR_indirect)

            PISM.verbPrintf(1, grid.com, "\nTest Jacobian Design Transpose (%d,%d):\n" % (i, j))

            view(r_global, r_viewer)
            view(JStarR, JStarR_viewer)
            view(JStarR_indirect, JStarR_indirect_viewer)
            view(d_JStarR, d_JStarR_viewer)

            PISM.logging.pause()


def test_linearization_transpose(ssarun):
    grid = ssarun.grid

    S = 250
    r_viewer = PETSc.Viewer().createDraw(title="r", size=S)
    TStarR_viewer = PETSc.Viewer().createDraw(title="TStarR", size=S)
    TStarR_indirect_viewer = PETSc.Viewer().createDraw(title="TStarR (ind)", size=S)
    d_TStarR_viewer = PETSc.Viewer().createDraw(title="d_TStarR_fd", size=S)

    ssarun.solver.linearize_at(zeta1)
    u = PISM.Vector1(grid, "")
    u.copy_from(ssarun.solver.solution())

    Td = PISM.Vector(grid, "")

    TStarR = PISM.Scalar(grid, "")

    TStarR_indirect = PISM.Scalar(grid, "")

    for (i, j) in grid.points():

        for k in range(2):

            r = PISM.Vector1(grid, "")
            r.set(0)
            with PISM.vec.Access(comm=r):
                if k == 0:
                    r[i, j].u = 1
                else:
                    r[i, j].v = 1

            ssarun.solver.apply_linearization_transpose(r, TStarR)

            r_global = PISM.Vector(grid, "")
            r_global.copy_from(r)

            for (k, l) in grid.points():
                with PISM.vec.Access(nocomm=TStarR_indirect):
                    d = PISM.Scalar1(grid, "")
                    d.set(0)
                    with PISM.vec.Access(comm=d):
                        d[k, l] = 1

                    ssarun.solver.apply_linearization(d, Td)

                    TStarR_indirect[k, l] = Td.get_vec().dot(r_global.get_vec())

            d_TStarR = PISM.Scalar(grid, "")

            d_TStarR.copy_from(TStarR)
            d_TStarR.add(-1, TStarR_indirect)

            PISM.verbPrintf(1, grid.com, "\nTest Linearization Transpose (%d,%d):\n" % (i, j))

            view(r_global, r_viewer)
            view(TStarR, TStarR_viewer)
            view(TStarR_indirect, TStarR_indirect_viewer)
            view(d_TStarR, d_TStarR_viewer)

            PISM.logging.pause()


# Main code starts here
if __name__ == "__main__":
    context = PISM.Context()
    config = context.config
    com = context.com
    PISM.set_abort_on_sigint(True)

    append_mode = False

    input_filename = config.get_string("input.file")
    inv_data_filename = PISM.OptionString("-inv_data", "inverse data file", input_filename).value()
    use_design_prior = config.get_flag("inverse.use_design_prior")
    design_var = PISM.OptionKeyword("-inv_ssa",
                                    "design variable for inversion",
                                    "tauc,hardav", "tauc").value()
    using_zeta_fixed_mask = config.get_flag("inverse.use_zeta_fixed_mask")

    ssarun = PISM.invert.ssa.SSAForwardRunFromInputFile(input_filename, inv_data_filename, design_var)
    ssarun.setup()

    vecs = ssarun.modeldata.vecs
    grid = ssarun.grid

    # Determine the prior guess for tauc/hardav. This can be one of
    # a) tauc/hardav from the input file (default)
    # b) tauc/hardav_prior from the inv_datafile if -inv_use_design_prior is set
    design_prior = createDesignVec(grid, design_var, '%s_prior' % design_var)
    long_name = design_prior.metadata().get_string("long_name")
    units = design_prior.metadata().get_string("units")
    design_prior.metadata().set_string("long_name",
                                       "best prior estimate for %s (used for inversion)" % long_name)
    if PISM.util.fileHasVariable(inv_data_filename, "%s_prior" % design_var) and use_design_prior:
        PISM.logging.logMessage("  Reading '%s_prior' from inverse data file %s.\n" % (design_var, inv_data_filename))
        design_prior.regrid(inv_data_filename, critical=True)
    else:
        if not PISM.util.fileHasVariable(input_filename, design_var):
            PISM.verbPrintf(1, com, "Initial guess for design variable is not available as '%s' in %s.\nYou can provide an initial guess in the inverse data file.\n" % (
                design_var, input_filename))
            exit(1)
        PISM.logging.logMessage("Reading '%s_prior' from '%s' in input file.\n" % (design_var, design_var))
        design = createDesignVec(grid, design_var)
        design.regrid(input_filename, True)
        design_prior.copy_from(design)

    if using_zeta_fixed_mask:
        if PISM.util.fileHasVariable(inv_data_filename, "zeta_fixed_mask"):
            zeta_fixed_mask = PISM.model.createZetaFixedMaskVec(grid)
            zeta_fixed_mask.regrid(inv_data_filename)
            vecs.add(zeta_fixed_mask)
        else:
            if design_var == 'tauc':
                logMessage(
                    "  Computing 'zeta_fixed_mask' (i.e. locations where design variable '%s' has a fixed value).\n" % design_var)
                zeta_fixed_mask = PISM.model.createZetaFixedMaskVec(grid)
                zeta_fixed_mask.set(1)
                mask = vecs.mask
                with PISM.vec.Access(comm=zeta_fixed_mask, nocomm=mask):
                    for (i, j) in grid.points():
                        if mask.grounded_ice(i, j):
                            zeta_fixed_mask[i, j] = 0
                vecs.add(zeta_fixed_mask)

                adjustTauc(vecs.mask, design_prior)
            elif design_var == 'hardav':
                pass
            else:
                raise NotImplementedError("Unable to build 'zeta_fixed_mask' for design variable %s.", design_var)

    # Convert design_prior -> zeta_prior
    zeta1 = PISM.Scalar2(grid, "")
    ssarun.designVariableParameterization().convertFromDesignVariable(design_prior, zeta1)

    ssarun.solver.linearize_at(zeta1)

    test_type = PISM.OptionKeyword("-inv_test", "",
                                   "j_design,j_design_transpose,lin,lin_transpose",
                                   "j_design")

    if not test_type.is_set():
        PISM.verbPrintf(1, com, "Must specify a test type via -inv_test\n")
        exit(1)
    else:
        test_type = test_type.value()

    if test_type == "j_design":
        test_j_design(ssarun)
    elif test_type == "j_design_transpose":
        test_j_design_transpose(ssarun)
    elif test_type == "lin":
        test_lin(ssarun)
    elif test_type == "lin_transpose":
        test_linearization_transpose(ssarun)
