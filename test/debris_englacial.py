#!/usr/bin/env python3
"""Verification of the englacial debris bookkeeping (`pism::debris::EnglacialTransport`):
advection, burial and melt-out in a single column, compared with exact solutions.

Emergence in the ablation zone: ice moves up at speed `a` while surface melt lowers the
surface at the same rate, so the thickness is steady. A uniform initial concentration
`C0` is advected toward the surface; the debris contained in the melted ice is released.
The exact cumulative melt-out after time `T` is `C0 a T` and the remaining englacial mass
is `C0 (H - a T)` (as long as the debris-free ice entering at the bed has not reached the
surface). The concentration front moving up from the bed is smeared by the scheme but
does not affect these integrals.

Burial in the accumulation zone: with a prescribed input rate `F` and accumulation at the
rate `b`, the buried mass after `T` is `F T rho_s`, and the thickness of the debris-loaded
layer at the top is `b T`.
"""

import unittest
import numpy as np
import PISM
import PISM.testing

ctx = PISM.Context()
ctx.log.set_threshold(1)

rho_s = 0.7 * 2650.0


def create_grid(Mz, Lz):
    params = PISM.GridParameters(ctx.config, 3, 3, 1.0, 1.0)
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


def uniform_3d(grid, value):
    result = PISM.Array3D(grid, "v", PISM.WITH_GHOSTS, grid.z(), 1)
    result.set(value)
    return result


def scalar(grid, name, value):
    result = PISM.Scalar(grid, name)
    result.set(value)
    return result


def column_mass(transport):
    "Mass per unit area of the column at (0, 0)."
    total = PISM.Scalar(transport.mass().grid(), "total")
    transport.column_mass(total)
    return PISM.testing.sample(total)


class MeltOut(unittest.TestCase):

    def test_emergence_and_melt_out(self):
        H = 80.0
        Mz = 21
        grid = create_grid(Mz, 100.0)
        geometry = create_geometry(grid, H)

        a = 1e-7                       # emergence velocity, m/s (about 3 m/yr)
        dt = 1e7                       # Courant number 0.2 for dz = 5 m
        N = 20                         # a N dt = 20 m of ice melted

        C0 = 2.0
        C = PISM.Array3D(grid, "C", PISM.WITHOUT_GHOSTS, grid.z())
        C.set(C0)

        scheme = PISM.TransportScheme3D.create(grid, "upwind", 1, False)
        transport = PISM.EnglacialTransport(grid, scheme, rho_s)
        transport.set_concentration(C, geometry.ice_thickness)

        mass_0 = column_mass(transport)
        np.testing.assert_almost_equal(mass_0, C0 * H)

        top_smb = scalar(grid, "top_smb", -a * dt)
        bottom_smb = scalar(grid, "bottom_smb", 0.0)
        rate = scalar(grid, "rate", 0.0)
        u, v, w = uniform_3d(grid, 0.0), uniform_3d(grid, 0.0), uniform_3d(grid, a)

        melted = 0.0
        for _ in range(N):
            transport.step(dt, geometry.ice_thickness, geometry.cell_type,
                           geometry.ice_thickness, geometry.cell_type,
                           top_smb, bottom_smb, rate, u, v, w)
            melted += PISM.testing.sample(transport.melt_out())
            self.assertEqual(PISM.testing.sample(transport.burial()), 0.0)
            self.assertEqual(PISM.testing.sample(transport.basal_loss()), 0.0)
            self.assertEqual(PISM.testing.sample(transport.ice_free_loss()), 0.0)

        T = N * dt
        # the far tail of the smeared front reaches the surface at the 1e-9 level
        np.testing.assert_allclose(melted, C0 * a * T, rtol=1e-6)
        np.testing.assert_allclose(column_mass(transport), C0 * (H - a * T), rtol=1e-6)

        # the debris-free ice that entered at the bed occupies [0, a T] (smeared by the
        # scheme); the concentration well above the front is unchanged
        result = PISM.Array3D(grid, "C", PISM.WITHOUT_GHOSTS, grid.z())
        transport.concentration(geometry.ice_thickness, result)
        with PISM.vec.Access(nocomm=result):
            profile = np.array(result.get_column(0, 0)) if grid.ctx().rank() == 0 else None
        if profile is not None:
            z = np.array(grid.z())
            self.assertLess(profile[0], 0.1 * C0)
            np.testing.assert_allclose(profile[(z > 3 * a * T) & (z < H - 5.0)], C0, rtol=1e-3)

    def test_basal_melt(self):
        "Basal melt removes debris at the bed; the rest keeps its height above the bed."
        H = 80.0
        grid = create_grid(21, 100.0)
        geometry = create_geometry(grid, H)

        C0 = 3.0
        C = PISM.Array3D(grid, "C", PISM.WITHOUT_GHOSTS, grid.z())
        C.set(C0)

        transport = PISM.EnglacialTransport(grid, PISM.TransportScheme3D.create(grid, "upwind", 1, False), rho_s)
        transport.set_concentration(C, geometry.ice_thickness)

        dH = 4.0
        thinner = create_geometry(grid, H - dH)
        zero = uniform_3d(grid, 0.0)

        transport.step(1.0, geometry.ice_thickness, geometry.cell_type,
                       thinner.ice_thickness, thinner.cell_type,
                       scalar(grid, "top", 0.0), scalar(grid, "bottom", -dH), scalar(grid, "rate", 0.0),
                       zero, zero, zero)

        np.testing.assert_allclose(PISM.testing.sample(transport.basal_loss()), C0 * dH, rtol=1e-10)
        np.testing.assert_allclose(column_mass(transport), C0 * (H - dH), rtol=1e-10)


class Burial(unittest.TestCase):

    def test_burial(self):
        # The initial and final surfaces sit on interfaces of the control volumes (2.5 m
        # and 52.5 m for 5 m spacing), so that melting the added ice returns exactly the
        # buried mass; in general the piecewise-constant representation mixes the buried
        # debris with the rest of the volume the surface is in.
        H0 = 47.5
        grid = create_grid(21, 100.0)

        b = 5e-8                       # accumulation rate, m/s: 0.5 m per step
        F = 1e-9                       # debris input rate, m of solid debris per second
        dt = 1e7
        N = 10                         # 5 m of ice added: one whole volume

        transport = PISM.EnglacialTransport(grid, PISM.TransportScheme3D.create(grid, "upwind", 1, False), rho_s)
        zero = uniform_3d(grid, 0.0)
        rate = scalar(grid, "rate", F)
        bottom = scalar(grid, "bottom", 0.0)
        top = scalar(grid, "top", b * dt)

        H = H0
        buried = 0.0
        for _ in range(N):
            old = create_geometry(grid, H)
            H += b * dt
            new = create_geometry(grid, H)
            transport.step(dt, old.ice_thickness, old.cell_type, new.ice_thickness, new.cell_type,
                           top, bottom, rate, zero, zero, zero)
            buried += PISM.testing.sample(transport.burial())
            self.assertEqual(PISM.testing.sample(transport.melt_out()), 0.0)

        T = N * dt
        np.testing.assert_allclose(buried, F * T * rho_s, rtol=1e-12)
        np.testing.assert_allclose(column_mass(transport), F * T * rho_s, rtol=1e-12)

        # all of it sits in the ice added since the start: melting that ice releases it
        # (note: keep the Geometry objects alive while their fields are in use)
        thicker = create_geometry(grid, H)
        thinner = create_geometry(grid, H0)
        transport.step(dt, thicker.ice_thickness, thicker.cell_type,
                       thinner.ice_thickness, thinner.cell_type,
                       scalar(grid, "top", H0 - H), bottom, scalar(grid, "rate", 0.0), zero, zero, zero)
        np.testing.assert_allclose(PISM.testing.sample(transport.melt_out()), F * T * rho_s, rtol=1e-10)
        np.testing.assert_allclose(column_mass(transport), 0.0, atol=1e-9)

    def test_ice_free(self):
        "A column that loses all its ice gives up its debris."
        grid = create_grid(11, 100.0)
        transport = PISM.EnglacialTransport(grid, PISM.TransportScheme3D.create(grid, "upwind", 1, False), rho_s)
        C = PISM.Array3D(grid, "C", PISM.WITHOUT_GHOSTS, grid.z())
        C.set(1.0)
        old = create_geometry(grid, 40.0)
        transport.set_concentration(C, old.ice_thickness)
        gone = create_geometry(grid, 0.0)
        zero = uniform_3d(grid, 0.0)
        transport.step(1.0, old.ice_thickness, old.cell_type, gone.ice_thickness, gone.cell_type,
                       scalar(grid, "top", -40.0), scalar(grid, "bottom", 0.0), scalar(grid, "rate", 0.0),
                       zero, zero, zero)
        np.testing.assert_allclose(PISM.testing.sample(transport.ice_free_loss()), 40.0, rtol=1e-12)
        self.assertEqual(column_mass(transport), 0.0)


if __name__ == "__main__":
    unittest.main()
