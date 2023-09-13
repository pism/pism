#!/usr/bin/env python3
"""This script runs verification tests of the Blatter stress balance solver that use
manufactured solutions from

Tezaur et al (2015) Albany/FELIX: a parallel, scalable and robust, finite element,
first-order Stokes approximation ice sheet solver built for advanced analysis
(doi:10.5194/gmd-8-1197-2015)

"""

from unittest import TestCase

import PISM
import PISM.util

import numpy as np

# Keep petsc4py from suppressing error messages
PISM.PETSc.Sys.popErrorHandler()

ctx = PISM.Context()
config = ctx.config

# do not ignore thin ice
config.set_number("geometry.ice_free_thickness_standard", 0.0)
config.set_number("stress_balance.ice_free_thickness_standard", 0.0)

# Set flow law parameters
config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")
config.set_number("stress_balance.blatter.Glen_exponent", 3.0)

# Set constants:
config.set_number("constants.ice.density", 910.0)
config.set_number("constants.standard_gravity", 9.81)

config_clean = PISM.DefaultConfig(ctx.com, "pism_config", "-config", ctx.unit_system)
config_clean.init_with_default(ctx.log)
config_clean.import_from(config)

def expt(ns, errors):
    "Compute the convergence rate using a polynomial fit."
    return -np.polyfit(np.log(ns), np.log(errors), 1)[0]

class TestXY(TestCase):
    """2D (x-y) verification test using a manufactured solution.

    u = exp(x) * sin(2 * pi * y)
    v = exp(x) * cos(2 * pi * y)

    on [0, 1] * [0, 1] with Dirichlet BC at x = {0,1}, y = {0,1}.

    Flat bed, no basal drag, constant ice thickness, isothermal Glen flow law with n == 3.

    See section 4.1 in Tezaur et al.

    The source term needed for the chosen manufactured solution is computed (and the code
    is generated) using SymPy.

    """
    def setUp(self):
        "Set PETSc options"
        opt = PISM.PETSc.Options()

        self.opt = opt

        self.A_old = config.get_number("flow_law.isothermal_Glen.ice_softness")
        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-4, "Pa-3 year-1", "Pa-3 s-1"))

        self.opts = {"-bp_snes_monitor_ratio": "",
                     "-bp_ksp_type": "preonly",
                     "-bp_pc_type": "lu",
                     }

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

    def tearDown(self):
        "Clear PETSc options"
        config.import_from(config_clean)

        for k, v in self.opts.items():
            self.opt.delValue(k)

    def inputs(self, N):
        "Allocate stress balance inputs for a given grid size"

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

        grid = PISM.Grid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.Scalar(grid, "tauc")
        yield_stress.set(0.0)

        geometry.bed_elevation.set(0.0)
        geometry.ice_thickness.set(1.0)
        geometry.sea_level_elevation.set(0.0)
        geometry.ensure_consistency(0.0)

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid):
        "Returns an array with the exact solution"
        exact = PISM.Vector(grid, "exact")

        with PISM.vec.Access(exact):
            for (i, j) in grid.points():
                x = grid.x(i)
                y = grid.y(j)
                exact[i, j] = PISM.blatter_xy_exact(x, y)

        return exact

    def error_norm(self, N):
        "Return the infinity norm of errors for the u and v components (as a tuple)."
        geometry, enthalpy, yield_stress = self.inputs(N)

        grid = enthalpy.grid()

        u_exact = self.exact_solution(grid)

        # no variation in the Z direction, so it's OK to use 2 vertical levels
        blatter_Mz = 2
        coarsening_factor = 1

        model = PISM.BlatterTestXY(grid, blatter_Mz, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        # compute the error
        error = PISM.Vector(grid, "error")
        error.copy_from(u_exact)
        error.add(-1.0, model.velocity())

        return error.norm(PISM.PETSc.NormType.NORM_INFINITY)

    def test(self):
        "Test that the convergence rate for the XY test is at least quadratic"
        Ns = [11, 21]
        norms = [self.error_norm(N) for N in Ns]

        norms_u = [n[0] for n in norms]
        norms_v = [n[1] for n in norms]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms_u)
        expt_v = expt(Ns, norms_v)

        print("U component conv. rate: dx^{}".format(expt_u))
        print("V component conv. rate: dx^{}".format(expt_v))

        # The convergence rate should be close to quadratic.
        assert expt_u >= 2.0
        assert expt_v >= 2.0

    def plot(self):
        Ns = [11, 21, 41, 81]
        try:
            self.setUp()
            norms = [self.error_norm(N) for N in Ns]
        finally:
            self.tearDown()

        # the domain is [0, 1]*[0, 1]*[0, 1]
        dxs = 1.0 / (np.array(Ns) - 1)

        norms_u = [x[0] for x in norms]
        norms_v = [x[1] for x in norms]

        p_u = np.polyfit(np.log(dxs), np.log(norms_u), 1)
        fit_u = np.exp(np.polyval(p_u, np.log(dxs)))

        p_v = np.polyfit(np.log(dxs), np.log(norms_v), 1)
        fit_v = np.exp(np.polyval(p_v, np.log(dxs)))

        from bokeh.plotting import figure, show, output_file
        from bokeh.layouts import gridplot

        output_file("test-xy.html")
        f = figure(title = "Blatter-Pattyn stress balance: verification test XY, x-component",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx", y_axis_label="max error")
        f.scatter(dxs, norms_u)
        f.line(dxs, norms_u, legend_label="max error", line_width=2)
        f.line(dxs, fit_u, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p_u[0]))

        g = figure(title = "Blatter-Pattyn stress balance: verification test XY, y-component",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx", y_axis_label="max error")
        g.scatter(dxs, norms_v)
        g.line(dxs, norms_v, legend_label="max error", line_width=2)
        g.line(dxs, fit_v, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p_v[0]))

        gp = gridplot([[f, g]])

        show(gp)

class TestXZ(TestCase):
    """2D (x-z) verification test using a manufactured solution.

    u = SIA + SSA exact solutions
    v = 0

    With sliding and a constant drag (beta), constant ice thickness, parabolic top surface
    and bed elevations, isothermal Glen flow law with n == 3.

    Domain: x in [-50km, 50km], ice thickness of 1km. The domain extent in the Y direction
    is chosen to set dy = dx with My == 3. This avoids possible element aspect ratio issues.

    See section 4.2 in Tezaur et al.

    - Uses Dirichlet BC at x == +-50km (unlike Tezaur et al who use a stress BC derived
      from the manufactured solution).

    - Uses basal BC which is a sum of the sliding BC and a correction computed using the
      manufactured solution.

    - Uses a stress BC at the top surface derived from the manufactured solution.

    All source terms and BC formulas are from Tezaur et al. (The C code is generated using
    SymPy.)

    This verification test uses the code implementing the basal boundary condition and so
    allows us to check its correctness (at least when beta is a constant).

    """
    def setUp(self):
        "Set PETSc options"

        self.opt = PISM.PETSc.Options()

        self.opts = {"-bp_snes_monitor_ratio": "",
                     "-bp_pc_type": "mg",
                     "-bp_pc_mg_levels": "1", # added here to make tearDown clean it up
                     "-bp_snes_ksp_ew": "",
                     }

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

        # the magnitude of the regularization constant affects the accuracy near the
        # "dome"
        config.set_number("flow_law.Schoof_regularizing_velocity", 1e-5)

        # Set sliding law parameters to make "tauc" equivalent to "beta"
        config.set_flag("basal_resistance.pseudo_plastic.enabled", True)
        config.set_number("basal_resistance.pseudo_plastic.q", 1.0)
        config.set_number("basal_resistance.pseudo_plastic.u_threshold",
                          PISM.util.convert(1.0, "m / s", "m / year"))

        # set ice softness
        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-16, "Pa-3 year-1", "Pa-3 s-1"))

        self.A = config.get_number("flow_law.isothermal_Glen.ice_softness")

        self.s0    = 2000.0     # m
        self.alpha = 4e-8       # 1/m
        self.H     = 1000.0     # m
        self.beta  = PISM.util.convert(1.0, "kPa year m-1", "Pa s m-1")

    def tearDown(self):
        "Clear PETSc options and configuration parameters"

        config.import_from(config_clean)

        for k, v in self.opts.items():
            self.opt.delValue(k)

    def inputs(self, N):
        P = PISM.GridParameters(config)

        # Domain: [-50e3, 50e3] * [-dx, dx] * [0, 1000]
        Lx = 50e3
        Mx = N
        dx = (2 * Lx) / (Mx - 1)

        P.Lx = 50e3
        P.Mx = int(N)
        P.x0 = 0.0

        P.Ly = dx
        P.My = 3
        P.y0 = 0.0

        # this vertical grid is used to store ice enthalpy, not ice velocity, so 2 levels
        # is enough
        P.z = PISM.DoubleVector([0.0, self.H])
        P.registration = PISM.CELL_CORNER
        P.periodicity = PISM.Y_PERIODIC
        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.Grid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.Scalar(grid, "tauc")

        yield_stress.set(self.beta)

        with PISM.vec.Access(geometry.bed_elevation):
            for (i, j) in grid.points():
                x = grid.x(i)
                geometry.bed_elevation[i, j] = self.s0 - self.H - self.alpha * x**2

        geometry.ice_thickness.set(self.H)

        # ensure that all ice is grounded
        geometry.sea_level_elevation.copy_from(geometry.bed_elevation)
        geometry.sea_level_elevation.shift(-1.0)

        geometry.ensure_consistency(0.0)

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid, bed, Z):
        "Returns an array with the exact solution"
        exact = PISM.Array3D(grid, "exact", PISM.WITHOUT_GHOSTS, Z)
        exact.metadata(0).long_name("x-component of the exact solution").units("m / s").output_units("m / year")

        rho = config.get_number("constants.ice.density")
        g = config.get_number("constants.standard_gravity")

        u = np.zeros_like(Z)
        with PISM.vec.Access([bed, exact]):
            for (i, j) in grid.points():
                x = grid.x(i)
                for k, z_sigma in enumerate(Z):
                    z = bed[i, j] + self.H * z_sigma
                    u[k] = PISM.blatter_xz_exact(x, z, self.A, rho, g,
                                                 self.s0, self.alpha, self.H, self.beta).u
                exact.set_column(i, j, u)

        return exact

    def error_norm(self, N, n_mg):
        "Return the infinity norm of errors for the u component."
        geometry, enthalpy, yield_stress = self.inputs(N)

        # set the number of multigrid levels
        self.opt.setValue("-bp_pc_mg_levels", n_mg)

        grid = enthalpy.grid()

        blatter_Mz = N
        coarsening_factor = 4

        model = PISM.BlatterTestXZ(grid, blatter_Mz, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        u_model_z = model.velocity_u_sigma().get_levels()

        u_model = PISM.Array3D(grid, "u_model", PISM.WITHOUT_GHOSTS, u_model_z)
        u_model.metadata(0).long_name("modeled velocity").units("m / s").output_units("m / year")
        u_model.copy_from(model.velocity_u_sigma())

        u_exact = self.exact_solution(grid, geometry.bed_elevation, u_model_z)

        # compute the error
        u_error = PISM.Array3D(grid, "error", PISM.WITHOUT_GHOSTS, u_model_z)
        u_error.copy_from(u_exact)
        u_error.add(-1.0, model.velocity_u_sigma())

        return u_error.norm(PISM.PETSc.NormType.NORM_INFINITY)[0]

    def test(self):
        "Test that the convergence rate for the XZ test is at least quadratic"

        # refinement path
        Ns = [11, 21]
        # number of MG levels to use for a particular grid size, assuming the coarsening
        # factor of 4
        mg_levels = [1, 2]

        norms = [self.error_norm(N, n_mg) for (N, n_mg) in zip(Ns, mg_levels)]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms)

        print("U component conv. rate: dx^{}".format(expt_u))

        # The convergence rate should be close to quadratic.
        assert expt_u >= 2.0

    def plot(self):
        Ns = [11, 21, 49, 129]
        mg_levels = [1, 2, 3, 4]

        try:
            self.setUp()
            norms = [self.error_norm(N, n_mg) for (N, n_mg) in zip(Ns, mg_levels)]
        finally:
            self.tearDown()

        # the domain is 100km long and 1km high
        dxs = 100e3 / (np.array(Ns) - 1)

        p = np.polyfit(np.log(dxs), np.log(norms), 1)
        fit = np.exp(np.polyval(p, np.log(dxs)))

        from bokeh.plotting import figure, show, output_file

        output_file("test-xz.html")
        f = figure(title = "Blatter-Pattyn stress balance: verification test XZ",
                   x_axis_type="log", y_axis_type="log",
                   x_axis_label="dx", y_axis_label="max error")
        f.scatter(dxs, norms)
        f.line(dxs, norms, legend_label="max error", line_width=2)
        f.line(dxs, fit, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p[0]))


        show(f)

class TestCFBC(TestCase):
    """Constant viscosity 2D (x-z) verification test checking the implementation of CFBC.

    on a square 0 <= x <= 1, -1 <= z <= 0 with sea level at z = 0 (fully submerged).

    Dirichet BC at x = 0, periodic in the Y direction, lateral BC at x = 1. No basal drag.

    """
    def setUp(self):
        "Set PETSc options"

        self.H = 1e3
        self.L = 1.0

        self.opt = PISM.PETSc.Options()

        self.opts = {"-bp_snes_monitor_ratio": "",
                     "-bp_snes_ksp_ew": ""}

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

        # constant viscocity
        n = 1.0
        config.set_number("stress_balance.blatter.Glen_exponent", n)

        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-3, "Pa-3 year-1", "Pa-3 s-1"))
        A = config.get_number("flow_law.isothermal_Glen.ice_softness")
        self.B = 1.0 / A              # note that n = 1

    def tearDown(self):
        "Clear PETSc options and configuration parameters"

        config.import_from(config_clean)

        for k, v in self.opts.items():
            self.opt.delValue(k)

    def inputs(self, N):
        P = PISM.GridParameters(config)

        # Domain: [0, 1] * [-dx, dx] * [-1, 0]
        Lx = 0.5 * self.L
        Mx = N
        dx = (2 * Lx) / (Mx - 1)

        P.Lx = Lx
        P.Mx = int(N)
        P.x0 = Lx

        P.Ly = dx
        P.My = 3
        P.y0 = 0.0

        # this vertical grid is used to store ice enthalpy, not ice velocity, so 2 levels
        # is enough
        P.z = PISM.DoubleVector([0.0, self.H])
        P.registration = PISM.CELL_CORNER
        P.periodicity = PISM.Y_PERIODIC
        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.Grid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.Scalar(grid, "tauc")
        # this value is not important: we use a compensatory term at the base instead of
        # the sliding law
        yield_stress.set(0.0)

        geometry.bed_elevation.set(-self.H)
        geometry.ice_thickness.set(self.H)
        geometry.ice_surface_elevation.set(0.0)
        geometry.cell_type.set(PISM.MASK_FLOATING)
        geometry.sea_level_elevation.set(0.0)

        # do *not* call geometry.ensure_consistency(): we want to keep surface elevation
        # at the sea levels

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid, bed, Z):
        "Returns an array with the exact solution"
        exact = PISM.Array3D(grid, "exact", PISM.WITHOUT_GHOSTS, Z)
        exact.metadata(0).long_name("x-component of the exact solution").units("m / s").output_units("m / year")

        rho_i = config.get_number("constants.ice.density")
        rho_w = config.get_number("constants.sea_water.density")
        g     = config.get_number("constants.standard_gravity")

        u = np.zeros_like(Z)
        with PISM.vec.Access([bed, exact]):
            for (i, j) in grid.points():
                x = grid.x(i)
                for k, z_sigma in enumerate(Z):
                    z = bed[i, j] + self.H * z_sigma
                    u[k] = PISM.blatter_xz_cfbc_exact(x, z, self.B, self.L, rho_i, rho_w, g).u
                exact.set_column(i, j, u)

        return exact

    def error_norm(self, N):
        "Return the infinity norm of errors for the u component."
        geometry, enthalpy, yield_stress = self.inputs(N)

        grid = enthalpy.grid()

        blatter_Mz = N
        coarsening_factor = 1

        model = PISM.BlatterTestCFBC(grid, blatter_Mz, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        u_model_z = model.velocity_u_sigma().get_levels()

        u_model = PISM.Array3D(grid, "u_model", PISM.WITHOUT_GHOSTS, u_model_z)
        u_model.metadata(0).long_name("modeled velocity").units("m / s").output_units("m / year")
        u_model.copy_from(model.velocity_u_sigma())

        u_exact = self.exact_solution(grid, geometry.bed_elevation, u_model_z)

        # compute the error
        u_error = PISM.Array3D(grid, "error", PISM.WITHOUT_GHOSTS, u_model_z)
        u_error.metadata(0).long_name("error").units("m / s").output_units("m / year")
        u_error.copy_from(u_exact)
        u_error.add(-1.0, model.velocity_u_sigma())

        return u_error.norm(PISM.PETSc.NormType.NORM_INFINITY)[0]

    def test(self):
        "Test that the convergence rate for the XZ-CFBC test is at least quadratic"

        # refinement path
        Ns = [11, 21]

        norms = [self.error_norm(N) for N in Ns]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms)

        print("U component conv. rate: dx^{}".format(expt_u))

        # The convergence rate should be close to quadratic.
        assert expt_u >= 2.0

    def plot(self):
        Ns = [11, 21, 41, 81]

        try:
            self.setUp()
            norms = [self.error_norm(N) for N in Ns]
        finally:
            self.tearDown()

        dxs = self.L / (np.array(Ns) - 1)

        p = np.polyfit(np.log(dxs), np.log(norms), 1)
        fit = np.exp(np.polyval(p, np.log(dxs)))

        from bokeh.plotting import figure, show, output_file

        output_file("test-xz-cfbc.html")
        f = figure(title = "Blatter-Pattyn stress balance: verification test XZ-CFBC",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx",
                   y_axis_label="max error")
        f.scatter(dxs, norms)
        f.line(dxs, norms, legend_label="max error", line_width=2)
        f.line(dxs, fit, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p[0]))


        show(f)

class TestXZvanderVeen(TestCase):
    def setUp(self):
        self.opt = PISM.PETSc.Options()

        self.opts = {"-bp_snes_monitor_ratio": "",
                     "-bp_ksp_type": "preonly",
                     "-bp_pc_type": "lu"}

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

        # Set sliding law parameters to make "tauc" equivalent to "beta"
        config.set_flag("basal_resistance.pseudo_plastic.enabled", True)
        config.set_number("basal_resistance.pseudo_plastic.q", 1.0)
        config.set_number("basal_resistance.pseudo_plastic.u_threshold",
                          PISM.util.convert(1.0, "m / s", "m / year"))

    def tearDown(self):
        for k in self.opts.keys():
            self.opt.delValue(k)
        config.import_from(config_clean)

    def create_grid(self, Mx=201):
        Lx = 1e5

        # compute dx and set Ly so that dy == dx
        dx = (2 * Lx) / (Mx - 1)
        dy = dx

        P = PISM.GridParameters(config)
        P.horizontal_size_from_options()
        P.Mx = Mx
        P.My = 3
        P.periodicity = PISM.Y_PERIODIC
        P.Lx = Lx
        P.x0 = 1.1 * Lx
        P.Ly = dy
        P.y0 = 0
        P.registration = PISM.CELL_CORNER
        P.z = PISM.DoubleVector([0, 1000])
        P.ownership_ranges_from_options(ctx.com.size)

        grid = PISM.Grid(ctx.ctx, P)

        return grid

    def compute(self, grid):
        tauc = PISM.Scalar(grid, "tauc")
        tauc.set(0.0)

        enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        enthalpy.set(0.0)

        Mz = 2
        coarsening_factor = 1
        model = PISM.BlatterTestvanderVeen(grid, Mz, coarsening_factor)

        geometry = PISM.Geometry(grid)

        # low enough to make it grounded
        sea_level = -100.0

        with PISM.vec.Access([geometry.ice_thickness, geometry.bed_elevation, tauc]):
            for (i, j) in grid.points():
                X = grid.x(i)
                geometry.ice_thickness[i, j] = model.H_exact(X)
                geometry.bed_elevation[i, j] = model.b_exact(X)
                tauc[i, j] = model.beta_exact(X)

        geometry.sea_level_elevation.set(sea_level)
        geometry.ensure_consistency(0.0)

        model.init()

        inputs = PISM.StressBalanceInputs()
        inputs.geometry = geometry
        inputs.basal_yield_stress = tauc
        inputs.enthalpy = enthalpy

        model.update(inputs, True)

        return model

    def error_norm(self, Mx):
        grid = self.create_grid(Mx)

        model = self.compute(grid)

        u_model = model.velocity_u_sigma().numpy()[1, :, :].mean(axis=1)
        u_exact = [model.u_exact(t).u for t in grid.x()]

        return np.max(np.fabs(u_model - u_exact))

    def test(self):
        "Check that the van der Veen flow line case converges"
        C = 10
        m = 2
        N = 7

        Mxs    = [C * 2**k + 1 for k in range(m, N)]
        errors = [self.error_norm(Mx) for Mx in Mxs]

        assert expt(Mxs, errors) >= 2.0

    def plot(self):
        C = 10
        m = 2
        N = 7
        Mxs = [C * 2**k + 1 for k in range(m, N)]

        try:
            self.setUp()
            errors = [self.error_norm(Mx) for Mx in Mxs]
        finally:
            self.tearDown()

        p = np.polyfit(np.log(Mxs), np.log(errors), 1)
        fit = np.exp(np.polyval(p, np.log(Mxs)))

        from bokeh.plotting import figure, show, output_file
        output_file("test-xz-van-der-Veen.html")
        f = figure(title="BP stress balance: verification test XZ (van der Veen profile)",
                   x_axis_type="log", y_axis_type="log",
                   x_axis_label="Mx", y_axis_label="max error")
        f.scatter(Mxs, errors)
        f.line(Mxs, fit, legend_label=f"O(dx^{-p[0]:1.2f})")
        show(f)

class TestXZHalfar(TestCase):
    def setUp(self):
        self.R0 = 750e3
        self.H0 = 3600.0

        self.opt = PISM.PETSc.Options()

        self.opts = {"-bp_snes_monitor_ratio": "",
                     "-bp_ksp_type": "cg",
                     "-bp_pc_type": "mg",
                     "-bp_pc_mg_levels": "1",
                     "-bp_snes_ksp_ew": ""}

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

    def tearDown(self):
        for k in self.opts.keys():
            self.opt.delValue(k)
        config.import_from(config_clean)

    def grid_whole(self, Mx):
        padding = 1e3
        Lx = self.R0 / 2 - padding
        x0 = self.R0 / 2

        return self.grid(Mx, Lx, x0)

    def grid_center(self, Mx):
        Lx = self.R0 / 4
        x0 = self.R0 / 2

        return self.grid(Mx, Lx, x0)

    def grid(self, Mx, Lx, x0):

        dx = (2 * Lx) / (Mx - 1)

        P = PISM.GridParameters(config)

        P.Lx = Lx
        P.Mx = Mx
        P.x0 = Lx

        P.Ly = dx
        P.My = 3
        P.y0 = 0.0

        # this vertical grid is used to store ice enthalpy, not ice velocity, so 2 levels
        # is enough
        P.z = PISM.DoubleVector([0.0, self.H0])
        P.registration = PISM.CELL_CORNER
        P.periodicity = PISM.Y_PERIODIC
        P.ownership_ranges_from_options(ctx.size)

        return PISM.Grid(ctx.ctx, P)

    def compute(self, grid, Mz, coarsening_factor):
        tauc = PISM.Scalar(grid, "tauc")
        tauc.set(0.0)

        enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        enthalpy.set(0.0)

        model = PISM.BlatterTestHalfar(grid, Mz, coarsening_factor)

        geometry = PISM.Geometry(grid)

        with PISM.vec.Access(geometry.ice_thickness):
            for (i, j) in grid.points():
                geometry.ice_thickness[i, j] = model.H_exact(grid.x(i))

        geometry.bed_elevation.set(0.0)
        geometry.sea_level_elevation.set(0)
        geometry.ensure_consistency(0.0)

        model.init()

        inputs = PISM.StressBalanceInputs()
        inputs.geometry = geometry
        inputs.basal_yield_stress = tauc
        inputs.enthalpy = enthalpy

        model.update(inputs, True)

        return model

    def exact_solution(self, model):
        "Returns an array with the exact solution"

        Z = model.velocity_u_sigma().get_levels()
        grid = model.grid()

        exact = PISM.Array3D(grid, "exact", PISM.WITHOUT_GHOSTS, Z)
        exact.metadata(0).long_name("x-component of the exact solution").units("m / s").output_units("m / year")

        u = np.zeros_like(Z)
        with PISM.vec.Access(exact):
            for (i, j) in grid.points():
                x = grid.x(i)
                H = model.H_exact(x)
                for k, z_sigma in enumerate(Z):
                    z = H * z_sigma
                    u[k] = model.u_exact(x, z)
                exact.set_column(i, j, u)

        return exact

    def error_norm(self, grid, Mz, mg_levels, coarsening_factor):
        self.opt.setValue("-bp_pc_mg_levels", str(mg_levels))
        model = self.compute(grid, Mz, coarsening_factor)

        u_exact = self.exact_solution(model)

        # compute the error
        Z = model.velocity_u_sigma().get_levels()
        u_error = PISM.Array3D(grid, "error", PISM.WITHOUT_GHOSTS, Z)
        u_error.copy_from(u_exact)
        u_error.add(-1.0, model.velocity_u_sigma())

        return u_error.norm(PISM.PETSc.NormType.NORM_INFINITY)[0]

    def test(self):
        "Check that the Halfar test case converges with z refinement"

        F = 2
        Mx = 51
        mg_levels = [5, 6]
        Mzs = [F**(m-1) + 1 for m in mg_levels]

        grid = self.grid_center(Mx)
        errors = [self.error_norm(grid, int(Mz), N, F) for Mz, N in zip(Mzs, mg_levels)]

        assert expt(Mzs, errors) >= 2.0

    def plot(self):
        self.plot_Mx()
        self.plot_Mz()

    def plot_Mx(self):
        # coarsening factor
        Mxs = [101, 201, 401]
        mg_levels = 3
        F = 8
        Mz = int(F**(mg_levels - 1)) + 1

        try:
            self.setUp()
            errors = [self.error_norm(self.grid_whole(Mx), Mz, mg_levels, F)
                      for Mx in Mxs]
        finally:
            self.tearDown()

        p = np.polyfit(np.log(Mxs), np.log(errors), 1)
        fit = np.exp(np.polyval(p, np.log(Mxs)))

        from bokeh.plotting import figure, show, output_file
        from bokeh.layouts import gridplot

        output_file("test-xz-Halfar-Mx.html")
        f = figure(title="BP stress balance: verification test XZ (Halfar dome)",
                   x_axis_type="log", y_axis_type="log",
                   x_axis_label="Mx", y_axis_label="max error")
        f.scatter(Mxs, errors)
        f.line(Mxs, fit, legend_label=f"O(Mx^{p[0]:1.2f})")

        show(f)

    def plot_Mz(self):
        # coarsening factor
        F = 2
        Mx = 51
        mg_levels = [5, 6, 7]
        Mzs = [int(F**(m-1)) + 1 for m in mg_levels]

        try:
            self.setUp()
            grid = self.grid_center(Mx)
            errors = [self.error_norm(grid, Mz, N, F)
                      for Mz, N in zip(Mzs, mg_levels)]
        finally:
            self.tearDown()

        p = np.polyfit(np.log(Mzs), np.log(errors), 1)
        fit = np.exp(np.polyval(p, np.log(Mzs)))

        from bokeh.plotting import figure, show, output_file
        output_file("test-xz-Halfar-Mz.html")
        f = figure(title="BP stress balance: verification test XZ (Halfar dome)",
                   x_axis_type="log", y_axis_type="log",
                   x_axis_label="Mz", y_axis_label="max error")
        f.scatter(Mzs, errors)
        f.line(Mzs, fit, legend_label=f"O(Mz^{p[0]:1.2f})")
        show(f)

if __name__ == "__main__":

    for test in [TestXY(),
                 TestXZ(),
                 TestCFBC(),
                 TestXZvanderVeen(),
                 TestXZHalfar()]:
        test.plot()
