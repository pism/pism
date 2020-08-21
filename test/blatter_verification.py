#!/usr/bin/env python3
"""This script runs verification tests of the Blatter stress balance solver that use
manufactured solutions.

"""
from unittest import TestCase

import PISM
import numpy as np
from numpy import exp, sin, cos, pi

ctx = PISM.Context()
config = ctx.config

class TestXY(TestCase):
    def exact(self, x, y):
        "Exact solution"
        u0 = exp(x) * sin(2 * pi * y)
        v0 = exp(x) * cos(2 * pi * y)
        return u0, v0

    def setUp(self):
        "Set PETSc options"
        opt = PISM.PETSc.Options()

        self.opt = opt

        # this verification test has no basal drag and the default preconditioner is
        # ineffective
        opt.setValue("-pc_type", "gamg")
        opt.setValue("-snes_monitor", "")
        opt.setValue("-ksp_monitor", "")

    def tearDown(self):
        "Clear PETSc options"
        opt = self.opt

        opt.delValue("-pc_type")
        opt.delValue("-snes_monitor")
        opt.delValue("-ksp_monitor")

    def inputs(self, N):
        "Allocate stress balance inputs for a given grid size"
        config.set_number("geometry.ice_free_thickness_standard", 0.0)
        config.set_number("stress_balance.ice_free_thickness_standard", 0.0)

        P = PISM.GridParameters(config)

        # Domain: [0,1] * [0, 1] * [0, 1]
        P.Lx = 0.5
        P.Ly = 0.5
        P.x0 = 0.5
        P.y0 = 0.5
        # This vertical grid is used for the enthalpy field only *and* we use an
        # isothermal flow law, so 2 levels is enough.
        P.z = PISM.DoubleVector([0.0, 1.0])
        P.registration = PISM.CELL_CORNER
        P.Mx = int(N)
        P.My = int(N)
        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.IceGrid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)

        geometry.bed_elevation.set(0.0)
        geometry.ice_thickness.set(1.0)
        geometry.sea_level_elevation.set(0.0)
        geometry.ensure_consistency(0.0)

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid):
        "Returns an array with the exact solution"
        exact = PISM.IceModelVec2V(grid, "exact", PISM.WITHOUT_GHOSTS)

        with PISM.vec.Access(exact):
            for (i, j) in grid.points():
                x = grid.x(i)
                y = grid.y(j)
                V = self.exact(x, y)
                exact[i, j].u = V[0]
                exact[i, j].v = V[1]

        return exact

    def error_norm(self, N):
        "Return the infinity norm of errors for the u and v components (as a tuple)."
        geometry, enthalpy, yield_stress = self.inputs(N)

        grid = enthalpy.grid()

        u_exact = self.exact_solution(grid)

        # no variation in the Z direction, so it's OK to use 2 vertical levels
        blatter_Mz = 2
        # no coarsening in the Z direction
        n_levels = 0
        coarsening_factor = 1

        model = PISM.BlatterTest1(grid, blatter_Mz, n_levels, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        # compute the error
        error = PISM.IceModelVec2V(grid, "error", PISM.WITHOUT_GHOSTS)
        error.copy_from(u_exact)
        error.add(-1.0, model.velocity())

        return error.norm_all(PISM.PETSc.NormType.NORM_INFINITY)

    def exact_xy_test(self):
        "Test that the convergence rate for the XY test is close to quadratic"
        Ns = [31, 61, 121]
        norms = [self.error_norm(N) for N in Ns]

        norms_u = [n[0] for n in norms]
        norms_v = [n[1] for n in norms]

        def expt(xs, ys):
            return -np.polyfit(np.log(xs), np.log(ys), 1)[0]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms_u)
        expt_v = expt(Ns, norms_v)

        print("U component conv. rate: dx^{}".format(expt_u))
        print("V component conv. rate: dx^{}".format(expt_v))

        # The convergence rate should be close to quadratic.
        assert np.fabs(expt_u - 2) < 0.1
        assert np.fabs(expt_v - 2) < 0.1
