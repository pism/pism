#!/usr/bin/env python3
"""Verification of the 3D transport schemes used for englacial debris
(`pism::TransportScheme3D`: first-order upwinding and MPDATA with and without the FCT
limiter) against an exact solution of the advection equation: translation of a Gaussian
by a uniform velocity in a slab of constant thickness, periodic in x and y, over whole
periods in the horizontal and a fraction of the slab thickness in the vertical (the blob
stays away from the bed and the surface, so the exact solution is the initial condition
shifted vertically).

The transported quantity is mass per unit area per control volume (see
`pism::debris::column`), so the exact solution is the concentration integrated over each
volume.
"""

import unittest
import numpy as np
import PISM

ctx = PISM.Context()
ctx.log.set_threshold(1)


def convergence_rate(dxs, errors):
    return np.polyfit(np.log(dxs), np.log(errors), 1)[0]


def interfaces(z):
    "Interfaces of the control volumes around the levels z (see column::interfaces)."
    zi = np.empty(len(z) + 1)
    zi[0] = 0.0
    zi[1:-1] = 0.5 * (z[:-1] + z[1:])
    zi[-1] = np.inf
    return zi


def create_grid(M, Mz, Lz):
    "Periodic [-1,1]^2 x [0, Lz] grid with equally spaced levels."
    params = PISM.GridParameters(ctx.config, M, M, 1.0, 1.0)
    params.periodicity = PISM.XY_PERIODIC
    params.z = PISM.DoubleVector(list(np.linspace(0.0, Lz, Mz)))
    params.ownership_ranges_from_options(ctx.config, ctx.size)
    return PISM.Grid(ctx.ctx, params)


def create_geometry(grid, H):
    geometry = PISM.Geometry(grid)
    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(-1.0)
    geometry.ice_thickness.set(H)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.ensure_consistency(0.0)
    return geometry


def gaussian_mass(grid, H, center, sigma):
    """Mass per unit area per control volume of a 3D Gaussian concentration, integrated
    analytically over each volume (the integral of exp(-z^2/2s^2) is an error function)."""
    from scipy.special import erf

    z = np.array(grid.z())
    zi = interfaces(z)
    x0, y0, z0 = center

    result = PISM.Array3D(grid, "mass", PISM.WITHOUT_GHOSTS, grid.z())

    with PISM.vec.Access(nocomm=result):
        for (i, j) in grid.points():
            x, y = grid.x(i), grid.y(j)
            horizontal = np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * sigma**2))
            column = []
            for k in range(len(z)):
                bottom = zi[k]
                top = min(zi[k + 1], H)
                if top <= bottom:
                    column.append(0.0)
                    continue
                s = sigma * np.sqrt(2.0)
                vertical = 0.5 * s * np.sqrt(np.pi) * (erf((top - z0) / s) - erf((bottom - z0) / s))
                column.append(horizontal * vertical)
            result.set_column(i, j, column)
    return result


def uniform_3d(grid, value):
    result = PISM.Array3D(grid, "v", PISM.WITH_GHOSTS, grid.z(), 1)
    result.set(value)
    return result


def total(x):
    "Total mass: sum over columns and levels times the cell area."
    column_sum = PISM.Scalar(x.grid(), "sum")
    PISM.sum_columns(x, 0.0, 1.0, column_sum)
    return PISM.sum(column_sum) * x.grid().cell_area()


def l1_error(x, exact):
    a = x.to_numpy()
    b = exact.to_numpy()
    if a is None:
        return None
    return np.sum(np.abs(a - b)) * x.grid().cell_area()


SCHEMES = {
    "upwind": ("upwind", 1, False, 0.4),
    "mpdata": ("mpdata", 2, False, 1.2),
    "mpdata-fct": ("mpdata", 2, True, 1.2),
}


class Translation(unittest.TestCase):

    def run_case(self, M, kind, N, fct):
        Lz = 2.0
        H = Lz                  # a slab filling the vertical grid
        Mz = M + 1              # dz = Lz / M = dx: the same resolution as horizontally

        grid = create_grid(M, Mz, Lz)
        geometry = create_geometry(grid, H)

        u, v, w = 1.0, 0.5, 0.1
        t_final = 4.0           # 2 periods in x, 1 in y; vertical shift w t = 0.4

        # the blob has to stay a few sigma away from the bed and the surface (closed faces)
        sigma = 0.2
        x = gaussian_mass(grid, H, (0.0, 0.0, 0.8), sigma)
        exact = gaussian_mass(grid, H, (0.0, 0.0, 0.8 + w * t_final), sigma)
        mass_0 = total(x)

        cfl_3d = PISM.max_timestep_cfl_3d(geometry.ice_thickness, geometry.cell_type, None,
                                          uniform_3d(grid, u), uniform_3d(grid, v),
                                          uniform_3d(grid, w))
        dz = Lz / (Mz - 1)
        dt_cfl = min(cfl_3d.dt_max.value(), dz / w)
        n_steps = int(np.ceil(t_final / (0.4 * dt_cfl)))
        dt = t_final / n_steps

        scheme = PISM.TransportScheme3D.create(grid, kind, N, fct)
        u3, v3, w3 = uniform_3d(grid, u), uniform_3d(grid, v), uniform_3d(grid, w)

        minimum = 0.0
        for _ in range(n_steps):
            scheme.update(dt, geometry.ice_thickness, geometry.cell_type, x, u3, v3, w3)
            x.copy_from(scheme.x())
        data = x.to_numpy()
        if data is not None:
            minimum = data.min()

        return l1_error(x, exact), (total(x) - mass_0) / mass_0, minimum

    def test_translation(self):
        for name, (kind, N, fct, min_rate) in SCHEMES.items():
            with self.subTest(scheme=name):
                Ms = [24, 48]
                errors, mass_errors, minima = zip(*[self.run_case(M, kind, N, fct) for M in Ms])

                rate = convergence_rate([2.0 / M for M in Ms], errors)
                self.assertGreater(rate, min_rate, f"{name}: errors {errors}, rate {rate:.2f}")
                for e in mass_errors:
                    self.assertLess(abs(e), 1e-12, f"{name}: mass not conserved ({e:g})")
                if fct or kind == "upwind":
                    for m in minima:
                        self.assertGreater(m, -1e-12, f"{name}: negative values ({m:g})")


class ClosedFaces(unittest.TestCase):
    """Mass cannot leave the ice: with velocities pointing out of the ice body nothing
    crosses the bed, the surface, or a face to an ice-free column."""

    def test_closed_domain(self):
        M, Mz, Lz = 12, 9, 4.0
        grid = create_grid(M, Mz, Lz)

        # a slab of thickness 2.5 in the middle third (in x) of the domain, ice-free elsewhere
        geometry = PISM.Geometry(grid)
        geometry.bed_elevation.set(0.0)
        geometry.sea_level_elevation.set(-1.0)
        geometry.ice_area_specific_volume.set(0.0)
        with PISM.vec.Access(nocomm=geometry.ice_thickness):
            for (i, j) in grid.points():
                geometry.ice_thickness[i, j] = 2.5 if abs(grid.x(i)) < 0.34 else 0.0
        geometry.ensure_consistency(0.0)

        x = PISM.Array3D(grid, "mass", PISM.WITHOUT_GHOSTS, grid.z())
        with PISM.vec.Access(nocomm=[x, geometry.ice_thickness]):
            for (i, j) in grid.points():
                x.set_column(i, j, [1.0] * Mz if geometry.ice_thickness[i, j] > 0 else [0.0] * Mz)
        mass_0 = total(x)

        for (u, v, w) in [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0),
                          (0.3, 0.2, 0.7)]:
            scheme = PISM.TransportScheme3D.create(grid, "mpdata", 2, True)
            y = x.duplicate()
            y.copy_from(x)
            for _ in range(20):
                scheme.update(0.01, geometry.ice_thickness, geometry.cell_type, y,
                              uniform_3d(grid, u), uniform_3d(grid, v), uniform_3d(grid, w))
                y.copy_from(scheme.x())
            self.assertLess(abs(total(y) - mass_0) / mass_0, 1e-12, f"velocity {(u, v, w)}")
            data = y.to_numpy()
            if data is not None:
                self.assertGreater(data.min(), -1e-12)


class Factory(unittest.TestCase):
    def test_unknown_scheme(self):
        grid = create_grid(4, 3, 1.0)
        with self.assertRaises(RuntimeError):
            PISM.TransportScheme3D.create(grid, "uno2", 1, False)
        with self.assertRaises(RuntimeError):
            PISM.MPDATA3(grid, 0, False)


if __name__ == "__main__":
    unittest.main()
