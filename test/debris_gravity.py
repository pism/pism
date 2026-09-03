#!/usr/bin/env python3
"""Verification of the gravitational transport of supraglacial debris
(`pism::debris::GravitationalTransport`) against the Barenblatt-Pattle solution.

With a flat ice surface the debris thickness obeys the porous-medium equation

    dh/dt = c div(h grad h) = (c/2) Laplacian(h^2),    c = mu (1 - phi) rho g,

whose radially symmetric similarity solution (Barenblatt, 1952; Pattle, 1959) in two
dimensions is

    h(r, tau) = tau^(-1/2) max(0, A - r^2 / (16 sqrt(tau))),    tau = c t / 2.

This is the same family of solutions as PISM's Halfar test for the SIA.
"""

import unittest
import numpy as np
import PISM

ctx = PISM.Context()
ctx.log.set_threshold(1)


def convergence_rate(dxs, errors):
    return np.polyfit(np.log(dxs), np.log(errors), 1)[0]


def barenblatt(r, tau, A):
    return tau**-0.5 * np.maximum(0.0, A - r**2 / (16.0 * np.sqrt(tau)))


def create_grid(M):
    return PISM.Grid.Shallow(ctx.ctx, 1.0, 1.0, 0.0, 0.0, M, M, PISM.CELL_CENTER, PISM.NOT_PERIODIC)


def create_geometry(grid):
    "Flat, ice-covered domain."
    geometry = PISM.Geometry(grid)
    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(-1.0)
    geometry.ice_thickness.set(100.0)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.ensure_consistency(0.0)
    return geometry


def set_barenblatt(h, tau, A):
    grid = h.grid()
    with PISM.vec.Access(nocomm=h):
        for (i, j) in grid.points():
            r = np.hypot(grid.x(i), grid.y(j))
            h[i, j] = barenblatt(r, tau, A)
    h.update_ghosts()


def l1_error(h, exact):
    diff = h.duplicate()
    diff.copy_from(h)
    diff.add(-1.0, exact)
    return diff.norm(PISM.PETSc.NormType.N1)[0] * h.grid().cell_area()


def total(h):
    return PISM.sum(h) * h.grid().cell_area()


class Barenblatt(unittest.TestCase):

    def run_case(self, M, cfl_ratio=0.25):
        grid = create_grid(M)
        geometry = create_geometry(grid)

        # c = mobility * solid_density * g = 1 (arbitrary units)
        g = ctx.config.get_number("constants.standard_gravity")
        solid_density = 1.0
        mobility = 1.0 / (solid_density * g)
        c = 1.0

        transport = PISM.GravitationalTransport(grid, mobility, solid_density, cfl_ratio)
        np.testing.assert_almost_equal(transport.coefficient(), c)

        A = 0.1
        tau_0, tau_1 = 0.01, 0.05   # support radius grows from 0.4 to 0.6

        h = PISM.Scalar1(grid, "debris_thickness")
        set_barenblatt(h, tau_0, A)
        exact = PISM.Scalar1(grid, "exact")
        set_barenblatt(exact, tau_1, A)
        mass_0 = total(h)

        dt = 2.0 * (tau_1 - tau_0) / c
        # several calls, as in a model run; the sub-cycling handles stability
        n_calls = 4
        substeps = 0
        for _ in range(n_calls):
            transport.step(dt / n_calls, geometry.cell_type, geometry.ice_surface_elevation, h)
            substeps += transport.substeps()

        return l1_error(h, exact), (total(h) - mass_0) / mass_0, PISM.min(h), substeps

    def test_convergence(self):
        Ms = [41, 81]
        errors, mass_errors, minima, substeps = zip(*[self.run_case(M) for M in Ms])

        rate = convergence_rate([2.0 / M for M in Ms], errors)
        # first order: the solution has a corner at the edge of its support
        self.assertGreater(rate, 0.7, f"errors {errors}, rate {rate:.2f}")
        for e in mass_errors:
            self.assertLess(abs(e), 1e-12)
        for m in minima:
            self.assertGreaterEqual(m, 0.0)
        for s in substeps:
            self.assertGreater(s, 4)

    def test_substep_independence(self):
        """Halving the stability ratio doubles the number of sub-steps; the result changes
        only by the (first-order) time discretization error."""
        e1, _, _, n1 = self.run_case(41, cfl_ratio=0.25)
        e2, _, _, n2 = self.run_case(41, cfl_ratio=0.125)
        self.assertGreater(n2, 1.8 * n1)
        self.assertLess(abs(e1 - e2), 0.3 * e1)


class ClosedMargin(unittest.TestCase):

    def test_no_flux_to_ice_free_cells(self):
        "Debris does not spread onto ice-free cells; mass is conserved on the ice."
        M = 31
        grid = create_grid(M)
        geometry = PISM.Geometry(grid)
        geometry.bed_elevation.set(0.0)
        geometry.sea_level_elevation.set(-1.0)
        geometry.ice_area_specific_volume.set(0.0)
        with PISM.vec.Access(nocomm=geometry.ice_thickness):
            for (i, j) in grid.points():
                geometry.ice_thickness[i, j] = 100.0 if grid.x(i) < 0.3 else 0.0
        geometry.ensure_consistency(0.0)

        g = ctx.config.get_number("constants.standard_gravity")
        transport = PISM.GravitationalTransport(grid, 1.0 / g, 1.0, 0.25)

        h = PISM.Scalar1(grid, "debris_thickness")
        set_barenblatt(h, 0.01, 0.1)
        with PISM.vec.Access(nocomm=[h, geometry.ice_thickness]):
            for (i, j) in grid.points():
                if geometry.ice_thickness[i, j] == 0.0:
                    h[i, j] = 0.0
        h.update_ghosts()
        mass_0 = total(h)

        for _ in range(5):
            transport.step(0.02, geometry.cell_type, geometry.ice_surface_elevation, h)

        self.assertLess(abs(total(h) - mass_0) / mass_0, 1e-12)
        data = h.to_numpy()
        if data is not None:
            x = np.array(grid.x())
            self.assertEqual(np.abs(data[:, x >= 0.3]).max(), 0.0)


if __name__ == "__main__":
    unittest.main()
