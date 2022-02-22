#!/usr/bin/env python3

from unittest import TestCase

import numpy as np
import scipy.integrate
import scipy.signal

import PISM

# silence models' initialization messages
PISM.Context().log.set_threshold(1)

class LingleClarkElastic(TestCase):
    @staticmethod
    def lrm(Mx, My, dx, dy):
        "Compute the load response matrix, taking advantage of its symmetry."

        ge = PISM.greens_elastic()

        def ge_integrand(eta, xi, dx, dy, p, q):
            xi_shift  = p * dx - xi
            eta_shift = q * dy - eta
            r         = np.sqrt(xi_shift * xi_shift + eta_shift * eta_shift)

            return ge(r)

        def f(dx, dy, p, q):
            return scipy.integrate.dblquad(ge_integrand,
                                           -dx/2.0, dx/2.0,
                                           lambda x: -dy / 2.0, lambda x: dy / 2.0,
                                           (dx, dy, p, q), epsrel=1e-8)[0]

        Mx2 = int(Mx) // 2
        My2 = int(My) // 2

        a = np.zeros((My, Mx))

        # top half
        for j in range(My2 + 1):
            # top left quarter
            for i in range(Mx2 + 1):
                p = Mx2 - i
                q = My2 - j

                a[j, i] = f(dx, dy, p, q)

            # top right quarter
            for i in range(Mx2 + 1, Mx):
                a[j, i] = a[j, 2 * Mx2 - i]

        # bottom half
        for j in range(My2 + 1, My):
            for i in range(Mx):
                a[j, i] = a[2 * My2 - j, i]

        return a

    @staticmethod
    def run_model(grid):
        geometry = PISM.Geometry(grid)

        bed_model = PISM.LingleClark(grid)

        bed_uplift = PISM.IceModelVec2S(grid, "uplift")

        # start with a flat bed, no ice, and no uplift
        geometry.bed_elevation.set(0.0)
        geometry.ice_thickness.set(0.0)
        geometry.sea_level_elevation.set(-1000.0) # everything is grounded
        geometry.ensure_consistency(0.0)

        bed_uplift.set(0.0)

        bed_model.bootstrap(geometry.bed_elevation, bed_uplift, geometry.ice_thickness,
                            geometry.sea_level_elevation)

        Mx2 = int(grid.Mx()) // 2
        My2 = int(grid.My()) // 2

        # add the load
        with PISM.vec.Access(nocomm=geometry.ice_thickness):
            for (i, j) in grid.points():
                # if i == Mx2 and j == My2:
                if abs(i - Mx2) < 2 and abs(j - My2) < 2:
                    geometry.ice_thickness[i, j] = 1000.0

        # dt of zero disables the viscous part of the model, so all we get is the elastic
        # response
        bed_model.step(geometry.ice_thickness, geometry.sea_level_elevation, 0)

        return (geometry.ice_thickness.numpy(),
                bed_model.total_displacement().numpy(),
                bed_model.elastic_load_response_matrix().numpy())

    def setUp(self):
        self.ctx = PISM.Context()

        self.elastic = self.ctx.config.get_flag("bed_deformation.lc.elastic_model")
        self.ctx.config.set_flag("bed_deformation.lc.elastic_model", True)

        self.size_factor = self.ctx.config.get_number("bed_deformation.lc.grid_size_factor")
        self.ctx.config.set_number("bed_deformation.lc.grid_size_factor", 2)

        Lx = 2000e3
        Mx = 11
        My = 2 * Mx             # non-square grid

        self.grid = PISM.IceGrid.Shallow(self.ctx.ctx, Lx, Lx,
                                         0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

        self.H, self.db_pism, self.lrm_pism = self.run_model(self.grid)

    def convolution_test(self):
        "Compare PISM's FFTW-based convolution to scipy.signal.fftconvolve()"
        rho = self.ctx.config.get_number("constants.ice.density")

        db_scipy = scipy.signal.fftconvolve(rho * self.H, self.lrm_pism, mode="same")

        np.testing.assert_allclose(self.db_pism, db_scipy, rtol=1e-12)

    def lrm_test(self):
        "Compare PISM's load response matrix to the one computed using scipy.integrate.dblquad()"

        dx = self.grid.dx()
        dy = self.grid.dy()

        Ny, Nx = self.lrm_pism.shape
        lrm_python = self.lrm(Nx, Ny, dx, dy)

        # This is a crappy relative tolerance. Oh well...
        np.testing.assert_allclose(self.lrm_pism, lrm_python, rtol=1e-2)

    def tearDown(self):
        # reset configuration parameters
        self.ctx.config.set_flag("bed_deformation.lc.elastic_model", self.elastic)
        self.ctx.config.set_number("bed_deformation.lc.grid_size_factor", self.size_factor)
