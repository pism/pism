#!/usr/bin/env python3
"""Verification of the 2D transport schemes used for supraglacial debris
(`pism::TransportScheme2D`: first-order upwinding, MPDATA with and without the FCT limiter,
UNO2, UNO3) against exact solutions of the advection equation:

- translation by a uniform velocity over one period of a periodic domain (the exact
  solution is the initial condition),
- solid-body rotation by one full revolution (same).

We check that the error decreases with resolution at (roughly) the expected rate, that mass
is conserved to rounding error, and that the limited schemes stay non-negative.
"""

import unittest
import numpy as np
import PISM

ctx = PISM.Context()
ctx.log.set_threshold(1)


def convergence_rate(dxs, errors):
    "Convergence rate p (error ~ dx^p) from a log-log fit."
    return np.polyfit(np.log(dxs), np.log(errors), 1)[0]


def create_grid(M, periodicity):
    "Unit square [-1,1]^2 with M x M cells."
    return PISM.Grid.Shallow(ctx.ctx, 1.0, 1.0, 0.0, 0.0, M, M, PISM.CELL_CENTER, periodicity)


def create_geometry(grid):
    "Ice everywhere, so that every cell is 'icy' and no face is closed."
    geometry = PISM.Geometry(grid)
    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(-1.0)
    geometry.ice_thickness.set(100.0)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.ensure_consistency(0.0)
    return geometry


def gaussian(grid, x0, y0, sigma):
    result = PISM.Scalar(grid, "x")
    with PISM.vec.Access(nocomm=result):
        for (i, j) in grid.points():
            x, y = grid.x(i), grid.y(j)
            result[i, j] = np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * sigma**2))
    return result


def uniform_velocity(grid, u, v):
    result = PISM.Vector(grid, "velocity")
    with PISM.vec.Access(nocomm=result):
        for (i, j) in grid.points():
            result[i, j] = PISM.Vector2d(u, v)
    result.update_ghosts()
    return result


def rotation_velocity(grid, omega):
    "Solid-body rotation around the origin: (u, v) = omega * (-y, x)."
    result = PISM.Vector(grid, "velocity")
    with PISM.vec.Access(nocomm=result):
        for (i, j) in grid.points():
            result[i, j] = PISM.Vector2d(-omega * grid.y(j), omega * grid.x(i))
    result.update_ghosts()
    return result


def total(x):
    return PISM.sum(x) * x.grid().cell_area()


def advect(scheme, geometry, x, velocity, t_final, cfl):
    "Advance x to t_final with a fixed fraction of the CFL time step. Returns (x, min)."
    # note: keep the CFLData object alive while reading its members
    cfl_data = PISM.max_timestep_cfl_2d(geometry.ice_thickness, geometry.cell_type, None,
                                        velocity)
    dt_cfl = cfl_data.dt_max.value()
    n_steps = int(np.ceil(t_final / (cfl * dt_cfl)))
    dt = t_final / n_steps

    minimum = PISM.min(x)
    for _ in range(n_steps):
        scheme.update(dt, geometry.cell_type, x, velocity)
        x.copy_from(scheme.x())
        minimum = min(minimum, PISM.min(x))

    return x, minimum


def l1_error(x, exact):
    diff = x.duplicate()
    diff.copy_from(x)
    diff.add(-1.0, exact)
    return diff.norm(PISM.PETSc.NormType.N1)[0] * x.grid().cell_area()


SCHEMES = {
    # name: (kind, N, nonoscillatory, minimum convergence rate, CFL number)
    #
    # UNO2/UNO3 use one-dimensional flux approximations in each direction, so for a flow
    # that is not aligned with the grid their accuracy degrades as the Courant number
    # grows (about first order at 0.5, second order at 0.25 for the diagonal translation
    # below); they are run with the smaller CFL number.
    "upwind": ("upwind", 1, False, 0.4, 0.5),
    "mpdata": ("mpdata", 2, False, 1.2, 0.5),
    "mpdata-fct": ("mpdata", 2, True, 1.2, 0.5),
    "uno2": ("uno2", 1, False, 0.9, 0.25),
    "uno3": ("uno3", 1, False, 0.9, 0.25),
}


class Translation(unittest.TestCase):
    """Uniform velocity, periodic domain, one period: exact = initial condition."""

    def run_case(self, M, kind, N, fct, cfl):
        grid = create_grid(M, PISM.XY_PERIODIC)
        geometry = create_geometry(grid)

        u, v = 1.0, 0.5
        # one period in x is 2 Lx / u; in y it is 2 Ly / v = 4: run for 4 time units so
        # that both directions complete whole periods (2 in x, 1 in y)
        t_final = 4.0

        x = gaussian(grid, 0.0, 0.0, 0.15)
        exact = x.duplicate()
        exact.copy_from(x)
        mass_0 = total(x)

        scheme = PISM.TransportScheme2D.create(grid, kind, N, fct)
        x, minimum = advect(scheme, geometry, x, uniform_velocity(grid, u, v), t_final, cfl)

        return l1_error(x, exact), (total(x) - mass_0) / mass_0, minimum

    def test_translation(self):
        for name, (kind, N, fct, min_rate, cfl) in SCHEMES.items():
            with self.subTest(scheme=name):
                Ms = [40, 80]
                errors, mass_errors, minima = zip(*[self.run_case(M, kind, N, fct, cfl) for M in Ms])

                rate = convergence_rate([2.0 / M for M in Ms], errors)
                self.assertGreater(rate, min_rate,
                                   f"{name}: errors {errors}, rate {rate:.2f}")
                for e in mass_errors:
                    self.assertLess(abs(e), 1e-12, f"{name}: mass not conserved ({e:g})")
                if fct or kind in ("upwind", "uno2", "uno3"):
                    for m in minima:
                        self.assertGreater(m, -1e-12, f"{name}: negative values ({m:g})")


class Rotation(unittest.TestCase):
    """Solid-body rotation by one revolution: exact = initial condition."""

    def run_case(self, M, kind, N, fct, cfl):
        grid = create_grid(M, PISM.NOT_PERIODIC)
        geometry = create_geometry(grid)

        omega = 2 * np.pi          # one revolution per time unit
        x = gaussian(grid, 0.4, 0.0, 0.12)
        exact = x.duplicate()
        exact.copy_from(x)
        mass_0 = total(x)

        scheme = PISM.TransportScheme2D.create(grid, kind, N, fct)
        x, minimum = advect(scheme, geometry, x, rotation_velocity(grid, omega), 1.0, cfl)

        return l1_error(x, exact), (total(x) - mass_0) / mass_0, minimum

    def test_rotation(self):
        for name, (kind, N, fct, min_rate, cfl) in SCHEMES.items():
            with self.subTest(scheme=name):
                Ms = [40, 80]
                errors, mass_errors, minima = zip(*[self.run_case(M, kind, N, fct, cfl) for M in Ms])

                rate = convergence_rate([2.0 / M for M in Ms], errors)
                self.assertGreater(rate, min_rate,
                                   f"{name}: errors {errors}, rate {rate:.2f}")
                # the blob never reaches the boundary, so mass is conserved
                for e in mass_errors:
                    self.assertLess(abs(e), 1e-10, f"{name}: mass not conserved ({e:g})")


class Factory(unittest.TestCase):
    def test_unknown_scheme(self):
        grid = create_grid(5, PISM.NOT_PERIODIC)
        with self.assertRaises(RuntimeError):
            PISM.TransportScheme2D.create(grid, "no-such-scheme", 1, False)
        with self.assertRaises(RuntimeError):
            PISM.TransportScheme2D.create(grid, "mpdata", 0, False)


if __name__ == "__main__":
    unittest.main()
