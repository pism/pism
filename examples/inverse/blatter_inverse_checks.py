#!/usr/bin/env python3
"""Consistency checks for the Blatter inverse forward problems.

Two checks are run for the selected design variable (``-design tauc`` or
``-design hardav``) on a small ISMIP-HOM C setup (see ``ismiphom_twin.py``):

1. **Adjoint (dot-product) test**: for random ``delta`` (design space) and
   ``w`` (surface-velocity space),

       <DF delta, w>  ==  <delta, DF^T w>

   where ``DF`` is ``apply_linearization`` and ``DF^T`` is
   ``apply_linearization_transpose``. This isolates the design Jacobian and
   its transpose from everything else. Agreement is limited by the linear
   solver tolerances, so the SNES/KSP tolerances are tightened here.

2. **Finite-difference gradient test**: for ``f(zeta) = 1/2 |F(zeta) - u_obs|^2``
   (plain nodal inner product) and ``g = DF^T (F(zeta) - u_obs)``,

       (f(zeta + h d) - f(zeta - h d)) / (2 h)  ==  <g, d>

   for a sequence of step sizes ``h``. The ratio of the two sides should
   approach 1.

Usage::

  python3 blatter_inverse_checks.py -design hardav [-Mx 21 -stress_balance.blatter.Mz 5 ...]
  mpiexec -n 4 python3 blatter_inverse_checks.py -design tauc

Both checks use the "exact" adjoint (KSPSolveTranspose with a direct
solver) unless ``-inverse.adjoint.method`` is given. The "approximate" and
"incomplete" (Picard) adjoints are approximations by construction, so with
them the identities above hold only approximately and the checks are
expected to report FAIL; use them only to gauge the size of that
approximation.
"""

import os
import sys

import numpy as np

import PISM
import PISM.invert.core
from PISM.util import convert

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ismiphom_twin as twin  # noqa: E402

seconds_per_year = convert(1.0, "year", "second")

forward_problems = {
    "tauc": PISM.IP_BlatterTaucForwardProblem,
    "hardav": PISM.IP_BlatterHardavForwardProblem,
}


def dot_scalar(grid, a, b):
    """Plain nodal inner product of two 2D scalar fields."""
    s = 0.0
    with PISM.vec.Access(nocomm=[a, b]):
        for (i, j) in grid.points():
            s += a[i, j] * b[i, j]
    return PISM.GlobalSum(grid.com, s)


def dot_vector(grid, a, b):
    """Plain nodal inner product of two 2D vector fields."""
    s = 0.0
    with PISM.vec.Access(nocomm=[a, b]):
        for (i, j) in grid.points():
            s += a[i, j].u * b[i, j].u + a[i, j].v * b[i, j].v
    return PISM.GlobalSum(grid.com, s)


def random_scalar(grid, name, seed, scale=1.0):
    """Ghostless scalar field with reproducible, rank-independent noise."""
    v = PISM.Scalar(grid, name)
    with PISM.vec.Access(comm=[v]):
        for (i, j) in grid.points():
            rng = np.random.default_rng(seed + 1000003 * i + 7919 * j)
            v[i, j] = scale * rng.standard_normal()
    return v


def random_vector(grid, name, seed, scale=1.0):
    v = PISM.Vector(grid, name)
    with PISM.vec.Access(comm=[v]):
        for (i, j) in grid.points():
            rng = np.random.default_rng(seed + 1000003 * i + 7919 * j)
            v[i, j].u = scale * rng.standard_normal()
            v[i, j].v = scale * rng.standard_normal()
    return v


def true_hardav(grid, L):
    """Smooth 'true' vertically-averaged hardness field."""
    config = grid.ctx().config()
    A = config.get_number("flow_law.isothermal_Glen.ice_softness")
    n = config.get_number("stress_balance.blatter.Glen_exponent")
    B0 = A ** (-1.0 / n)

    hardav = PISM.Scalar(grid, "hardav")
    hardav.metadata(0).long_name("vertically-averaged ice hardness")
    hardav.metadata(0).set_units_without_validation("Pa s^(1/n)")
    omega = 2.0 * np.pi / L
    with PISM.vec.Access(comm=[hardav]):
        for (i, j) in grid.points():
            x, y = grid.x(i), grid.y(j)
            hardav[i, j] = B0 * (1.0 + 0.3 * np.sin(omega * x) * np.sin(omega * y))
    return hardav


def misfit(grid, u, u_obs, residual):
    """f = 1/2 |u - u_obs|^2 (nodal); also stores u - u_obs in `residual`."""
    residual.copy_from(u)
    residual.add(-1.0, u_obs)
    return 0.5 * dot_vector(grid, residual, residual)


def main():
    ctx = PISM.Context()
    com = ctx.com
    config = ctx.config
    PISM.set_abort_on_sigint(True)

    design_var = PISM.OptionKeyword("-design", "design variable",
                                    "tauc,hardav", "tauc").value()
    L = 80e3

    twin.set_constants(config)
    config.set_string("inverse.design.param", "exp")
    # Tight forward tolerances: the checks compare quantities that are only
    # equal up to the linear/nonlinear solver accuracy.
    config.set_number("stress_balance.blatter.relative_convergence", 1e-12)

    opts = PISM.PETSc.Options()
    for key, value in [("-bp_snes_rtol", "1e-12"), ("-bp_snes_atol", "1e-14"),
                       ("-bp_snes_stol", "1e-14"), ("-bp_snes_max_it", "200"),
                       ("-bp_ksp_type", "preonly"), ("-bp_pc_type", "lu"),
                       ("-bp_snes_ksp_ew", "0"),
                       ("-inv_adj_ksp_type", "preonly"), ("-inv_adj_pc_type", "lu")]:
        if not opts.hasName(key):
            opts.setValue(key, value)
    if ctx.size > 1:
        for key in ("-bp_pc_factor_mat_solver_type", "-inv_adj_pc_factor_mat_solver_type"):
            if not opts.hasName(key):
                opts.setValue(key, "mumps")
    if not PISM.OptionString("-inverse.adjoint.method", "").is_set():
        config.set_string("inverse.adjoint.method", "exact")

    grid = twin.create_grid(config, ctx, "C", L)
    geometry, enthalpy, yield_stress = twin.init_geometry(grid, "C", L)
    enthalpy.set(0.0)

    Mz = int(config.get_number("stress_balance.blatter.Mz"))
    coarsening = int(config.get_number("stress_balance.blatter.coarsening_factor"))

    design_param = PISM.invert.core.createDesignVariableParam(config, design_var)

    # "True" design field, and a perturbed one to linearize about
    if design_var == "tauc":
        true_field = yield_stress
    else:
        true_field = true_hardav(grid, L)

    zeta_true = PISM.Scalar2(grid, "zeta_true")
    design_param.convertFromDesignVariable(true_field, zeta_true)

    solver = forward_problems[design_var](grid, Mz, coarsening, design_param)
    solver.init()

    adjoint_method = config.get_string("inverse.adjoint.method")
    PISM.verbPrintf(1, com, "\nBlatter inverse checks: design variable '%s', "
                    "adjoint method '%s', grid %dx%d, Mz %d\n"
                    % (design_var, adjoint_method, grid.Mx(), grid.My(), Mz))
    if adjoint_method != "exact":
        PISM.verbPrintf(1, com, "NOTE: the '%s' adjoint is an approximation; "
                        "the checks below are exact only for 'exact'.\n"
                        % adjoint_method)

    # Synthetic observations from the true field
    reason = solver.linearize_at(zeta_true)
    if reason.failed():
        raise RuntimeError("forward solve (true field) failed: %s" % reason.description())
    u_obs = PISM.Vector(grid, "u_obs")
    u_obs.copy_from(solver.solution())

    # Linearize about a perturbed design: zeta = zeta_true + 0.3 * noise
    zeta = PISM.Scalar2(grid, "zeta")
    zeta.copy_from(zeta_true)
    zeta.add(0.3, random_scalar(grid, "noise", seed=11))
    zeta.update_ghosts()

    reason = solver.linearize_at(zeta)
    if reason.failed():
        raise RuntimeError("forward solve failed: %s" % reason.description())

    # ------------------------------------------------------------------
    # 1. Adjoint (dot-product) test
    # ------------------------------------------------------------------
    delta = random_scalar(grid, "delta", seed=23)
    w = random_vector(grid, "w", seed=37, scale=1.0 / seconds_per_year)

    DF_delta = PISM.Vector(grid, "DF_delta")
    solver.apply_linearization(delta, DF_delta)

    DFt_w = PISM.Scalar(grid, "DFt_w")
    solver.apply_linearization_transpose(w, DFt_w)

    lhs = dot_vector(grid, DF_delta, w)
    rhs = dot_scalar(grid, delta, DFt_w)
    rel = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)
    PISM.verbPrintf(1, com,
                    "\nAdjoint test: <DF d, w> = %.12e  <d, DF^T w> = %.12e  "
                    "rel. difference = %.3e\n" % (lhs, rhs, rel))

    # ------------------------------------------------------------------
    # 2. Finite-difference gradient test
    # ------------------------------------------------------------------
    residual = PISM.Vector(grid, "residual")
    f0 = misfit(grid, solver.solution(), u_obs, residual)

    gradient = PISM.Scalar(grid, "gradient")
    solver.apply_linearization_transpose(residual, gradient)

    d = random_scalar(grid, "d", seed=41)
    g_dot_d = dot_scalar(grid, gradient, d)

    PISM.verbPrintf(1, com, "\nGradient test: f(zeta) = %.12e, <g, d> = %.12e\n"
                    % (f0, g_dot_d))
    PISM.verbPrintf(1, com, "  %10s %22s %22s %12s\n"
                    % ("h", "central difference", "<g, d>", "ratio"))

    zeta_h = PISM.Scalar2(grid, "zeta_h")
    ratios = []
    for h in [1e-1, 1e-2, 1e-3, 1e-4]:
        values = []
        for sign in (1.0, -1.0):
            zeta_h.copy_from(zeta)
            zeta_h.add(sign * h, d)
            zeta_h.update_ghosts()
            reason = solver.linearize_at(zeta_h)
            if reason.failed():
                raise RuntimeError("forward solve failed: %s" % reason.description())
            values.append(misfit(grid, solver.solution(), u_obs, residual))
        fd = (values[0] - values[1]) / (2.0 * h)
        ratio = fd / g_dot_d if g_dot_d != 0.0 else float("nan")
        ratios.append(ratio)
        PISM.verbPrintf(1, com, "  %10.1e %22.12e %22.12e %12.8f\n"
                        % (h, fd, g_dot_d, ratio))

    ok = rel < 1e-6 and min(abs(r - 1.0) for r in ratios) < 1e-3
    PISM.verbPrintf(1, com, "\nRESULT (%s): %s\n" % (design_var, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
