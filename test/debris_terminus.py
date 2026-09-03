#!/usr/bin/env python3
"""Tests of the removal of supraglacial debris at the glacier margin
(`pism::debris::TerminusRemoval`, equation 23 in Verhaegen and Huybrechts, 2026) and of
the debris mass budget of the complete transport model in a steady state where the input
from the terrain balances the removal at the margin.
"""

import os
import unittest
import numpy as np
import PISM
from PISM.testing import filename as tmp_name

ctx = PISM.Context()
ctx.log.set_threshold(1)
config = ctx.config

g = config.get_number("constants.standard_gravity")


def create_grid(Mx, My, Lx, Ly):
    "Cell-centered, non-periodic grid."
    return PISM.Grid.Shallow(ctx.ctx, Lx, Ly, 0.0, 0.0, Mx, My, PISM.CELL_CENTER, PISM.NOT_PERIODIC)


def strip_geometry(grid, x_margin, slope):
    """Ice for x < x_margin on a bed sloping down in x; flat in y. The first column of
    cells is an ice-free headwall (higher than the ice): PISM's ghosts wrap around the
    domain, so ice touching the domain boundary would "see" the far end of the domain as
    its neighbor, and an ice-free cell lower than the ice surface counts as foreland."""
    geometry = PISM.Geometry(grid)
    geometry.sea_level_elevation.set(-1e4)
    geometry.ice_area_specific_volume.set(0.0)
    with PISM.vec.Access(nocomm=[geometry.bed_elevation, geometry.ice_thickness]):
        for (i, j) in grid.points():
            x = grid.x(i)
            geometry.bed_elevation[i, j] = -slope * x if i > 0 else 1000.0
            geometry.ice_thickness[i, j] = 100.0 if (i > 0 and x < x_margin) else 0.0
    geometry.ensure_consistency(0.0)
    return geometry


class Removal(unittest.TestCase):

    def setUp(self):
        self.dx = 25.0
        Mx, My = 40, 5
        self.grid = create_grid(Mx, My, 0.5 * Mx * self.dx, 0.5 * My * self.dx)
        self.slope = 0.1
        # ice for i <= 29
        self.x_margin = self.grid.x(29) + 0.5 * self.dx
        self.geometry = strip_geometry(self.grid, self.x_margin, self.slope)
        self.rho_s = 0.7 * 2650.0
        self.mobility = 1e-4 / PISM.util.convert(1.0, "year", "second")

    def rate(self, gamma, h_d_value=1.0):
        removal = PISM.TerminusRemoval(self.grid, gamma, self.mobility, self.rho_s)

        h_d = PISM.Scalar2(self.grid, "debris_thickness")
        with PISM.vec.Access(nocomm=[h_d, self.geometry.ice_thickness]):
            for (i, j) in self.grid.points():
                h_d[i, j] = h_d_value if self.geometry.ice_thickness[i, j] > 0 else 0.0
        h_d.update_ghosts()

        dt = 1e6
        removal.step(dt, self.geometry.cell_type, self.geometry.ice_surface_elevation, h_d)

        rates = removal.removal_rate().to_numpy()
        thickness = h_d.to_numpy()
        return removal, rates, thickness, dt

    def test_margin_only(self):
        "Only the margin cell loses debris, at the rate of equation 23."
        removal, rates, thickness, dt = self.rate(gamma=self.dx)
        self.assertEqual(removal.upstream_cells(), 0)

        if rates is None:
            return

        # slope of the debris surface toward the ice-free neighbor: the surface drops by the
        # ice thickness plus the debris thickness plus the bed slope over one cell
        s = (100.0 + 1.0 + self.slope * self.dx) / self.dx
        expected = self.mobility * self.rho_s * g * 1.0 * s / self.dx

        self.assertEqual(np.count_nonzero(rates), self.grid.My())
        np.testing.assert_allclose(rates[:, 29], min(expected, 1.0 / dt), rtol=1e-12)
        np.testing.assert_allclose(thickness[:, 29], 1.0 - rates[0, 29] * dt, rtol=1e-12)
        np.testing.assert_allclose(thickness[:, 1:29], 1.0)

    def test_upstream_averaging(self):
        "With Gamma = 3 dx the thickness and slope are averaged over two more cells up-glacier."
        removal, rates, thickness, dt = self.rate(gamma=3 * self.dx)
        self.assertEqual(removal.upstream_cells(), 2)

        if rates is None:
            return

        s_margin = (100.0 + 1.0 + self.slope * self.dx) / self.dx
        s_up = self.slope      # the debris surface up-glacier slopes with the bed
        s_eff = (s_margin + 2 * s_up) / 3
        expected = self.mobility * self.rho_s * g * 1.0 * s_eff / (3 * self.dx)

        np.testing.assert_allclose(rates[:, 29], min(expected, 1.0 / dt), rtol=1e-12)

    def test_too_large_gamma(self):
        with self.assertRaises(RuntimeError):
            PISM.TerminusRemoval(self.grid, 4 * self.dx, self.mobility, self.rho_s)


class SteadyState(unittest.TestCase):
    """Debris added to the surface at one cell spreads downslope and leaves at the margin;
    in the steady state the output equals the input and the budget closes."""

    def setUp(self):
        self.dx = 25.0
        Mx, My = 40, 5
        self.grid = create_grid(Mx, My, 0.5 * Mx * self.dx, 0.5 * My * self.dx)
        self.geometry = strip_geometry(self.grid, self.grid.x(29) + 0.5 * self.dx, 0.1)

        # debris input at one cell in the middle of the strip
        self.input_file = tmp_name("debris_input")
        rate = PISM.Scalar(self.grid, "debris_input_rate")
        rate.metadata(0).long_name("debris input rate").units("m s^-1")
        self.F = 1e-9
        with PISM.vec.Access(nocomm=rate):
            for (i, j) in self.grid.points():
                rate[i, j] = self.F if (i == 20 and j == 2) else 0.0
        output = PISM.util.prepare_output(self.input_file)
        output.define_variable(rate.metadata())
        rate.write(output)
        output.close()

        config.set_string("debris.transport.input.file", self.input_file)
        config.set_string("debris.transport.supraglacial.scheme", "upwind")
        config.set_string("debris.transport.englacial.scheme", "upwind")
        # with the reference mobility the diffusive time scale over the strip is ~1e5
        # years; speed things up
        self.mobility = config.get_number("debris.transport.mobility")
        config.set_number("debris.transport.mobility", 0.1)

    def tearDown(self):
        os.remove(self.input_file)
        config.set_string("debris.transport.input.file", "")
        config.set_string("debris.transport.supraglacial.scheme", "mpdata")
        config.set_string("debris.transport.englacial.scheme", "mpdata")
        config.set_number("debris.transport.mobility", self.mobility)

    def test_budget(self):
        model = PISM.DebrisTransport(self.grid)
        model.init(self.geometry)

        zero3 = PISM.Array3D(self.grid, "zero", PISM.WITH_GHOSTS, self.grid.z(), 1)
        zero3.set(0.0)
        zero = PISM.Scalar(self.grid, "zero")
        zero.set(0.0)

        inputs = PISM.DebrisInputs()
        inputs.geometry = self.geometry
        inputs.u3 = zero3
        inputs.v3 = zero3
        inputs.w3 = zero3
        inputs.top_surface_mass_balance = zero
        inputs.bottom_surface_mass_balance = zero

        dt = PISM.util.convert(0.5, "year", "second")
        t = 0.0
        cell_area = self.grid.cell_area()
        input_flux = self.F * cell_area * model.solid_density()

        for n in range(4000):
            model.update(inputs, t, dt)
            t += dt
            # the budget closes at every step
            self.assertLess(abs(model.conservation_error()), 1e-9 * max(model.total_mass(), 1.0))
            np.testing.assert_allclose(model.input_last_step() / dt, input_flux, rtol=1e-12)

        # steady state: output balances input
        output_flux = model.output_last_step() / dt
        self.assertLess(abs(output_flux - input_flux) / input_flux, 1e-3)
        self.assertEqual(model.lost_last_step(), 0.0)
        self.assertEqual(model.englacial_mass(), 0.0)


if __name__ == "__main__":
    unittest.main()
