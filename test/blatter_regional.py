#!/usr/bin/env python3
"""Tests of the "regional mode" support in the Blatter stress balance solver.

These use small synthetic setups and drive the solver directly (no IceModel).

- Faces of elements at non-periodic domain edges are not treated as calving fronts in
  the regional mode.

- Velocity prescribed via bc_mask/bc_values (e.g. in the "no model" strip) is honored
  on all vertical levels.

- The driving stress in the "no model" strip uses the stored geometry
  (no_model_ice_thickness, no_model_surface_elevation): zero stored thickness means
  zero driving stress (regional.zero_gradient), and stored geometry equal to the current
  one reproduces the non-regional result.
"""

from unittest import TestCase

import PISM
import PISM.util

import numpy as np

# Keep petsc4py from suppressing error messages
PISM.PETSc.Sys.popErrorHandler()

ctx = PISM.Context()
config = ctx.config

config_clean = PISM.config_from_options(ctx.com, ctx.unit_system)
config_clean.import_from(config)

# number of grid points in each horizontal direction
N = 21
# number of vertical levels in the Blatter solver
MZ = 5

year = PISM.util.convert(1.0, "year", "seconds")


class Setup:
    """Allocates a grid, geometry and the other stress balance inputs for a synthetic
    slab and runs the Blatter solver with and without the regional mode."""

    def __init__(self, bed, thickness, sea_level, tauc=5e5, Lx=5e3, periodicity=PISM.NOT_PERIODIC):
        """bed and thickness are functions of (x, y) in meters.

        The domain is [0, 2*Lx] x [0, 2*Lx] with the cell-corner grid registration.
        """
        P = PISM.GridParameters(config, N, N, Lx, Lx)
        P.x0 = Lx
        P.y0 = Lx
        P.periodicity = periodicity
        # the vertical grid is used by the enthalpy field only; the flow law is isothermal
        P.z = PISM.DoubleVector([0.0, 5000.0])
        P.registration = PISM.CELL_CORNER
        P.ownership_ranges_from_options(ctx.config, ctx.size)

        self.grid = PISM.Grid(ctx.ctx, P)
        grid = self.grid

        self.geometry = PISM.Geometry(grid)

        with PISM.vec.Access([self.geometry.bed_elevation, self.geometry.ice_thickness]):
            for (i, j) in grid.points():
                x, y = grid.x(i), grid.y(j)
                self.geometry.bed_elevation[i, j] = bed(x, y)
                self.geometry.ice_thickness[i, j] = thickness(x, y)

        self.geometry.sea_level_elevation.set(sea_level)
        self.geometry.ensure_consistency(0.0)

        self.enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # the value is irrelevant: the flow law is isothermal
        self.enthalpy.set(1e5)

        self.yield_stress = PISM.Scalar(grid, "tauc")
        self.yield_stress.set(tauc)

        # regional inputs: by default the whole domain is "modeled"
        self.no_model_mask = PISM.Scalar2(grid, "no_model_mask")
        self.no_model_mask.set(0.0)

        self.no_model_thickness = PISM.Scalar(grid, "thkstore")
        self.no_model_thickness.copy_from(self.geometry.ice_thickness)

        self.no_model_surface = PISM.Scalar2(grid, "usurfstore")
        self.no_model_surface.copy_from(self.geometry.ice_surface_elevation)

        # Dirichlet BC: not used unless set_bc() is called
        self.bc_mask = None
        self.bc_values = None

    def set_no_model(self, in_strip, thickness=None, surface=None):
        """Mark grid points where in_strip(i, j) is True as part of the "no model" strip
        and set the stored geometry there (defaults: current geometry)."""
        with PISM.vec.Access([self.no_model_mask, self.no_model_thickness, self.no_model_surface,
                              self.geometry.ice_thickness, self.geometry.ice_surface_elevation]):
            for (i, j) in self.grid.points():
                if in_strip(i, j):
                    self.no_model_mask[i, j] = 1.0
                    if thickness is not None:
                        self.no_model_thickness[i, j] = thickness
                    if surface is not None:
                        self.no_model_surface[i, j] = surface

        self.no_model_mask.update_ghosts()
        self.no_model_surface.update_ghosts()

    def set_bc(self, in_strip, u, v):
        "Prescribe the velocity (u, v) (m/s) at grid points where in_strip(i, j) is True."
        self.bc_mask = PISM.Scalar(self.grid, "bc_mask")
        self.bc_mask.set(0.0)
        self.bc_values = PISM.Vector(self.grid, "_bc")
        self.bc_values.set(0.0)

        with PISM.vec.Access([self.bc_mask, self.bc_values]):
            for (i, j) in self.grid.points():
                if in_strip(i, j):
                    self.bc_mask[i, j] = 1.0
                    self.bc_values[i, j] = PISM.Vector2d(u, v)

    def solve(self, regional):
        "Run the Blatter solver; return (vertically averaged velocity, u_sigma, v_sigma)."
        coarsening_factor = 1
        model = PISM.Blatter(self.grid, MZ, coarsening_factor, regional)
        model.init()

        inputs = PISM.StressBalanceInputs()
        inputs.geometry = self.geometry
        inputs.basal_yield_stress = self.yield_stress
        inputs.enthalpy = self.enthalpy

        if regional:
            inputs.no_model_mask = self.no_model_mask
            inputs.no_model_ice_thickness = self.no_model_thickness
            inputs.no_model_surface_elevation = self.no_model_surface

        if self.bc_mask is not None:
            inputs.bc_mask = self.bc_mask
            inputs.bc_values = self.bc_values

        model.update(inputs, True)

        velocity = PISM.Vector(self.grid, "velocity")
        velocity.copy_from(model.velocity())

        return velocity, model.velocity_u_sigma(), model.velocity_v_sigma()


def max_norm(v):
    return v.norm(PISM.PETSc.NormType.NORM_INFINITY)


def difference(a, b):
    result = PISM.Vector(a.grid(), "difference")
    result.copy_from(a)
    result.add(-1.0, b)
    return result


class TestRegional(TestCase):

    def setUp(self):
        "Set PETSc options and flow law parameters"
        config.set_number("geometry.ice_free_thickness_standard", 0.0)
        config.set_number("stress_balance.ice_free_thickness_standard", 0.0)

        config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")
        config.set_number("stress_balance.blatter.Glen_exponent", 3.0)
        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-16, "Pa-3 year-1", "Pa-3 s-1"))
        config.set_flag("stress_balance.blatter.use_eta_transform", True)

        self.opts = {"-bp_snes_monitor_ratio": "",
                     "-bp_ksp_type": "preonly",
                     "-bp_pc_type": "lu",
                     }

        self.opt = PISM.PETSc.Options()
        for k, v in self.opts.items():
            self.opt.setValue(k, v)

    def tearDown(self):
        "Clear PETSc options"
        config.import_from(config_clean)

        for k, v in self.opts.items():
            self.opt.delValue(k)

    def test_domain_edges(self):
        """A uniform grounded slab with the bed below sea level: no driving stress, so
        the only forcing comes from the calving-front stress BC at domain edges.

        Without the regional mode this makes the ice move; in the regional mode domain
        edges are not calving fronts and the velocity is zero.
        """
        setup = Setup(bed=lambda x, y: -500.0,
                      thickness=lambda x, y: 1000.0,
                      sea_level=0.0)

        v_global, _, _ = setup.solve(regional=False)
        v_regional, _, _ = setup.solve(regional=True)

        u_global = max_norm(v_global)
        u_regional = max_norm(v_regional)
        print("max |u|: global = {} m/year, regional = {} m/year".format(u_global[0] * year,
                                                                          u_regional[0] * year))

        assert max(u_global) > 1e-10
        assert max(u_regional) < 1e-16

    def sloped_slab(self, periodicity=PISM.NOT_PERIODIC):
        "A slab with a sloped bed and uniform thickness; all ice is above sea level."
        return Setup(bed=lambda x, y: 1000.0 - 0.05 * x,
                     thickness=lambda x, y: 1000.0,
                     sea_level=-5000.0,
                     periodicity=periodicity)

    def test_periodic_grid(self):
        """In the regional mode the Blatter mesh is not periodic even if PISM's grid is
        (grid.periodicity is "xy" by default): the result on a periodic grid is the same
        as the non-regional result on a non-periodic grid.

        Without this, elements "wrapping around" the domain would see the jump in the
        bed elevation between the two edges of the sloped slab."""
        v_reference, _, _ = self.sloped_slab().solve(regional=False)
        assert max_norm(v_reference)[0] * year > 1.0

        v_periodic, _, _ = self.sloped_slab(PISM.XY_PERIODIC).solve(regional=True)

        d = max(max_norm(difference(v_periodic, v_reference)))
        print("max difference: {} m/year".format(d * year))
        assert d <= 1e-12 * max(max_norm(v_reference))

    @staticmethod
    def in_strip(i, j):
        "True in the three columns at the left and right edges of the domain"
        return i < 3 or i > N - 4

    def test_dirichlet_bc(self):
        "Velocity prescribed in the no-model strip is honored on all vertical levels."
        u0 = 10.0 / year
        v0 = -3.0 / year

        setup = self.sloped_slab()
        setup.set_no_model(self.in_strip)
        setup.set_bc(self.in_strip, u0, v0)

        velocity, u_sigma, v_sigma = setup.solve(regional=True)

        # 3D velocity components as arrays indexed by [j, i, k]
        u3 = u_sigma.to_numpy()
        v3 = v_sigma.to_numpy()

        grid = setup.grid
        interior = PISM.Vector(grid, "interior")
        interior.set(0.0)
        bc_error = 0.0
        with PISM.vec.Access([velocity, interior]):
            for (i, j) in grid.points():
                if self.in_strip(i, j):
                    U = velocity[i, j]
                    bc_error = max(bc_error, abs(U.u - u0), abs(U.v - v0))
                    bc_error = max(bc_error,
                                   np.max(np.abs(u3[j, i, :] - u0)),
                                   np.max(np.abs(v3[j, i, :] - v0)))
                else:
                    interior[i, j] = velocity[i, j]

        print("max BC error: {} m/year".format(bc_error * year))
        assert bc_error * year < 1e-10
        # ice in the modeled area moves down the slope
        assert max_norm(interior)[0] * year > 1.0

    def test_zero_gradient(self):
        """Zero stored thickness in the no-model strip (regional.zero_gradient) means zero
        driving stress there. With the whole domain in the strip nothing drives the flow,
        so the velocity is zero in spite of the sloped surface."""
        setup = self.sloped_slab()
        setup.set_no_model(lambda i, j: True, thickness=0.0, surface=0.0)

        v_regional, _, _ = setup.solve(regional=True)

        print("max |u| = {} m/year".format(max_norm(v_regional)[0] * year))
        assert max(max_norm(v_regional)) < 1e-16

    def test_stored_geometry(self):
        """Stored geometry equal to the current one reproduces the non-regional result,
        both with and without the no-model strip."""
        setup = self.sloped_slab()
        v_global, _, _ = setup.solve(regional=False)
        assert max_norm(v_global)[0] * year > 1.0

        v_regional, _, _ = setup.solve(regional=True)
        d = max(max_norm(difference(v_regional, v_global)))
        print("max difference (no strip): {} m/year".format(d * year))
        assert d <= 1e-12 * max(max_norm(v_global))

        setup.set_no_model(lambda i, j: True)
        v_regional, _, _ = setup.solve(regional=True)
        d = max(max_norm(difference(v_regional, v_global)))
        print("max difference (whole domain in the strip): {} m/year".format(d * year))
        assert d <= 1e-12 * max(max_norm(v_global))

    def test_zero_gradient_strip(self):
        """The zero-gradient treatment removes the driving stress in the no-model strip.

        Ice in the strip still moves because of the longitudinal coupling with the
        modeled area (the Blatter solver includes membrane stresses), but with the
        driving stress removed from all elements touching the strip the ice at the edge
        of the domain is much slower than with the stored geometry equal to the current
        one."""
        strip_width = 5

        def in_strip(i, j):
            return i < strip_width or i > N - 1 - strip_width

        def max_speed(velocity, condition):
            result = 0.0
            with PISM.vec.Access(velocity):
                for (i, j) in velocity.grid().points():
                    if condition(i, j):
                        U = velocity[i, j]
                        result = max(result, np.hypot(U.u, U.v))
            return result

        at_edge = lambda i, j: i < 2 or i > N - 3
        modeled = lambda i, j: not in_strip(i, j)

        # stored geometry equal to the current one
        setup = self.sloped_slab()
        setup.set_no_model(in_strip)
        v_stored, _, _ = setup.solve(regional=True)

        # zero stored thickness and surface elevation
        setup = self.sloped_slab()
        setup.set_no_model(in_strip, thickness=0.0, surface=0.0)
        v_zero, _, _ = setup.solve(regional=True)

        edge_stored = max_speed(v_stored, at_edge)
        edge_zero = max_speed(v_zero, at_edge)
        interior = max_speed(v_zero, modeled)

        print("max speed at the domain edge: stored geometry = {} m/year,"
              " zero gradient = {} m/year; modeled area = {} m/year".format(
                  edge_stored * year, edge_zero * year, interior * year))

        assert interior * year > 1.0
        assert edge_zero < 0.5 * edge_stored
