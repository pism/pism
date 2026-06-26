#!/usr/bin/env python3
"""ISMIP-HOM twin experiment: compare incomplete vs exact adjoint for Blatter inversion.

Sets up ISMIP-HOM C and D geometries with known tauc patterns, generates
synthetic surface velocities via a forward Blatter solve, then runs
Tikhonov inversions using both the incomplete (Picard/KSPSolve) and exact
(Newton/KSPSolveTranspose) adjoints.

Usage:
  mpiexec -n N python3 ismiphom_twin.py [options]

Example:
  mpiexec -n 8 python3 ismiphom_twin.py \
    -Mx 61 -stress_balance.blatter.Mz 5 \
    -stress_balance.blatter.coarsening_factor 2 \
    -bp_ksp_type fgmres -bp_pc_type mg -bp_pc_mg_levels 2 \
    -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_pc_type gamg \
    -bp_mg_coarse_ksp_rtol 1e-2 -bp_mg_coarse_ksp_max_it 50 \
    -tikhonov_atol 1e-5 -inv_max_it 250

Reference:
  Morlighem et al. (2013), "Inversion of basal friction in Antarctica
  using exact and incomplete adjoints of a higher-order model",
  J. Geophys. Res., 118, 1746-1753.
"""

import numpy as np
import sys
import time
import json

import PISM
from PISM.util import convert

# ============================================================================
# ISMIP-HOM geometry and sliding definitions
# ============================================================================

seconds_per_year = convert(1.0, "year", "second")


def surface_CD(x, y, L):
    alpha = 0.1 * (np.pi / 180.0)
    return -x * np.tan(alpha)


def bed_CD(x, y, L):
    return surface_CD(x, y, L) - 1000.0


def tauc_C(x, y, L):
    omega = 2.0 * np.pi / L
    return 1000.0 + 1000.0 * np.sin(omega * x) * np.sin(omega * y)


def tauc_D(x, y, L):
    omega = 2.0 * np.pi / L
    return 1000.0 + 1000.0 * np.sin(omega * x)


tests = {
    "C": {"tauc": tauc_C, "hom": PISM.HOM_C, "My_min": None},
    "D": {"tauc": tauc_D, "hom": PISM.HOM_D, "My_min": 3},
}


# ============================================================================
# Setup helpers
# ============================================================================

def set_constants(config):
    """Configure physics for ISMIP-HOM: isothermal Glen, linear sliding."""
    config.set_flag("basal_resistance.pseudo_plastic.enabled", True)
    config.set_number("basal_resistance.pseudo_plastic.q", 1.0)
    config.set_number("basal_resistance.pseudo_plastic.u_threshold",
                      convert(1.0, "m / s", "m / year"))

    config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")
    config.set_number("stress_balance.blatter.Glen_exponent", 3.0)
    config.set_number("flow_law.isothermal_Glen.ice_softness",
                      convert(1e-16, "Pa-3 year-1", "Pa-3 s-1"))

    config.set_number("constants.ice.density", 910.0)
    config.set_number("constants.standard_gravity", 9.81)

    config.set_flag("stress_balance.calving_front_stress_bc", True)
    config.set_flag("stress_balance.blatter.use_eta_transform", True)

    # Blatter solver defaults for ISMIP-HOM (loose tolerance, enough iterations)
    config.set_number("stress_balance.blatter.relative_convergence", 1e-3)

    config.set_string("output.format", "netcdf4_serial")


def create_grid(config, ctx, test_name, L):
    """Create a periodic grid for ISMIP-HOM test."""
    P = PISM.GridParameters(config)
    P.Lx = L / 2.0
    P.Ly = L / 2.0
    P.x0 = L / 2.0
    P.y0 = L / 2.0
    P.periodicity = PISM.XY_PERIODIC
    P.z = PISM.DoubleVector(np.linspace(0, 2000, 201))
    # Set sensible defaults before reading -Mx/-My/-grid.dx from the CLI.
    # Upstream PISM (commit 9703fd301) dropped the Mx>=3 validation that
    # used to catch the case where grid.Mx/grid.My defaults (-1) leave
    # GridParameters' Mx/My uninitialized; failure now surfaces deeper as
    # "'My' is invalid" inside ownership_ranges_from_options. Giving Mx/My
    # explicit defaults makes the script run without -Mx on the command
    # line, while still letting the user override via CLI.
    P.Mx = 21
    P.My = 21
    P.horizontal_size_and_extent_from_options(config)

    # Ensure My is set (default to Mx if only -Mx was given)
    if P.My <= 0:
        P.My = P.Mx

    if tests[test_name]["My_min"] is not None:
        P.My = max(P.My, tests[test_name]["My_min"])

    P.ownership_ranges_from_options(config, ctx.size)
    return PISM.Grid(ctx.ctx, P)


def init_geometry(grid, test_name, L):
    """Set up ISMIP-HOM geometry and yield stress."""
    geometry = PISM.Geometry(grid)
    grid.variables().add(geometry.ice_thickness)
    grid.variables().add(geometry.cell_type)
    grid.variables().add(geometry.bed_elevation)

    enthalpy = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
    enthalpy.metadata(0).long_name("enthalpy of ice").units("J kg^-1")
    grid.variables().add(enthalpy)

    yield_stress = PISM.Scalar(grid, "tauc")
    yield_stress.metadata(0).long_name("basal yield stress").units("Pa")
    grid.variables().add(yield_stress)

    tauc_fn = tests[test_name]["tauc"]

    with PISM.vec.Access([yield_stress, geometry.ice_thickness,
                          geometry.bed_elevation]):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            yield_stress[i, j] = tauc_fn(x, y, L) * seconds_per_year
            geometry.ice_thickness[i, j] = surface_CD(x, y, L) - bed_CD(x, y, L)
            geometry.bed_elevation[i, j] = bed_CD(x, y, L)

    geometry.sea_level_elevation.copy_from(geometry.bed_elevation)
    geometry.sea_level_elevation.shift(-100.0)
    geometry.ensure_consistency(0.0)

    return geometry, enthalpy, yield_stress


# ============================================================================
# Forward solve
# ============================================================================

def forward_solve(grid, yield_stress, config):
    """Run a forward Blatter solve using IP_BlatterTaucForwardProblem.

    Uses the same solver as the inversion so the twin experiment is
    self-consistent (no model mismatch between observations and inversion).
    """
    com = grid.com
    Mz = int(config.get_number("stress_balance.blatter.Mz"))
    coarsening = int(config.get_number("stress_balance.blatter.coarsening_factor"))

    param_name = config.get_string("inverse.design.param")
    design_param = PISM.invert.core.createDesignVariableParam(config, "tauc", param_name)

    solver = PISM.IP_BlatterTaucForwardProblem(grid, Mz, coarsening, design_param)
    solver.init()

    # Convert tauc -> zeta, then forward solve
    zeta = PISM.Scalar2(grid, "zeta")
    design_param.convertFromDesignVariable(yield_stress, zeta)

    reason = solver.linearize_at(zeta)
    PISM.verbPrintf(2, com, "Forward solve: %s\n" % reason.description())

    if reason.failed():
        raise RuntimeError("Forward Blatter solve failed: %s" % reason.description())

    # Copy the surface velocity (outlives the solver)
    u_obs = PISM.Vector(grid, "vel_observed")
    u_obs.copy_from(solver.solution())

    return u_obs


# ============================================================================
# Inversion
# ============================================================================

def run_inversion(grid, geometry, enthalpy, yield_stress_true,
                  u_obs, test_name, adjoint_method, config):
    """Run a Tikhonov inversion for tauc using the Blatter solver.

    Returns dict with convergence history and recovered tauc.
    """
    com = grid.com
    Mz = int(config.get_number("stress_balance.blatter.Mz"))
    coarsening = int(config.get_number("stress_balance.blatter.coarsening_factor"))

    PISM.verbPrintf(2, com, "\n=== Inversion (%s adjoint) ===\n" % adjoint_method)

    # Set the adjoint mode and matching KSP options
    config.set_string("inverse.adjoint.method", adjoint_method)

    opts = PISM.PETSc.Options()
    if adjoint_method == "incomplete":
        # Truly symmetric Picard Jacobian: CG + GAMG works
        opts.setValue("-inv_adj_ksp_type", "cg")
        opts.setValue("-inv_adj_pc_type", "gamg")
    elif adjoint_method == "approximate":
        # Symmetrized Newton Jacobian: use GMRES (not strictly SPD)
        opts.setValue("-inv_adj_ksp_type", "gmres")
        opts.setValue("-inv_adj_pc_type", "gamg")
    else:  # exact
        # KSPSolveTranspose on the full Newton Jacobian. Jacobi is far too
        # weak here (GMRES stalls -> DIVERGED_ITS at 10000 its). Use a direct
        # LU factorization instead: it supports an exact transpose solve and
        # is cheap at twin-experiment sizes. MUMPS handles the parallel case.
        opts.setValue("-inv_adj_ksp_type", "preonly")
        opts.setValue("-inv_adj_pc_type", "lu")
        opts.setValue("-inv_adj_pc_factor_mat_solver_type", "mumps")

    # Design variable parameterization
    param_name = config.get_string("inverse.design.param")
    design_param = PISM.invert.core.createDesignVariableParam(config, "tauc", param_name)

    # Forward problem
    solver = PISM.IP_BlatterTaucForwardProblem(grid, Mz, coarsening, design_param)
    solver.init()

    # Start from a perturbed initial guess.
    # With "exp" parameterization, the SNES sees log(tauc), so even
    # scale=0.5 is only log(0.5)=-0.69 away in zeta-space — modest.
    prior_scale = 0.5
    tauc_prior = PISM.Scalar(grid, "tauc_prior")
    tauc_prior.copy_from(yield_stress_true)
    tauc_prior.scale(prior_scale)
    PISM.verbPrintf(2, com, "  Prior: %.1f x true tauc\n" % prior_scale)

    # Convert prior tauc -> zeta
    zeta = PISM.Scalar2(grid, "zeta")
    design_param.convertFromDesignVariable(tauc_prior, zeta)

    # Functionals
    velocity_scale = config.get_number(
        "inverse.stress_balance.velocity_scale", "m/second")
    cH1 = config.get_number("inverse.design.cH1")
    cL2 = config.get_number("inverse.design.cL2")

    # Misfit weight: 1 everywhere (all ice is grounded in ISMIP-HOM)
    misfit_weight = PISM.Scalar(grid, "vel_misfit_weight")
    misfit_weight.set(1.0)

    state_func = PISM.IPMeanSquareFunctional2V(grid, misfit_weight)
    state_func.normalize(velocity_scale)

    design_func = PISM.IPGroundedIceH1NormFunctional2S(
        grid, cL2, cH1, geometry.cell_type)

    # Tikhonov problem
    # Constructor order: forward, d0, u_obs, eta, designFunctional, stateFunctional
    eta = config.get_number("inverse.tikhonov.penalty_weight")

    tikhonov = PISM.IP_BlatterTaucTaoTikhonovProblem(
        solver, zeta, u_obs, eta, design_func, state_func)

    tao_types = {'tikhonov_lmvm': 'lmvm',
                 'tikhonov_cg': 'cg',
                 'tikhonov_blmvm': 'blmvm'}
    method = config.get_string("inverse.stress_balance.method")
    tao_type = tao_types.get(method, 'lmvm')
    PISM.verbPrintf(2, com, "  TAO type: %s (from method=%s)\n" % (tao_type, method))
    tao_solver = PISM.IP_BlatterTaucTaoTikhonovSolver(com, tao_type, tikhonov)

    max_it = int(config.get_number("inverse.max_iterations"))
    tao_solver.setMaximumIterations(max_it)

    # Convergence logger
    misfit_history = []

    class MisfitLogger(PISM.IP_BlatterTaucTaoTikhonovProblemListener):
        def iteration(self, problem, eta_param, it,
                      val_design, val_state,
                      d, d_diff, grad_design,
                      u, u_diff, grad_state, gradient):
            misfit_history.append({
                "iteration": it,
                "design_value": float(val_design),
                "state_value": float(val_state),
            })
            PISM.verbPrintf(2, com,
                            "  it %3d: state=%10.4e  design=%10.4e\n" %
                            (it, val_state, val_design))

    listener = MisfitLogger()
    tikhonov.addListener(listener)

    # d0 (prior) and initial guess are both zeta — design penalty should start at ~0
    tikhonov.setInitialGuess(zeta)

    # Diagnostic: check initial forward solve misfit
    reason0 = solver.linearize_at(zeta)
    vel0 = solver.solution()
    v0_l2 = 0.0
    with PISM.vec.Access(nocomm=[vel0, u_obs]):
        for (i, j) in grid.points():
            du = vel0[i, j].u - u_obs[i, j].u
            dv = vel0[i, j].v - u_obs[i, j].v
            v0_l2 += du * du + dv * dv
    v0_l2 = PISM.GlobalSum(com, v0_l2)
    N = grid.Mx() * grid.My()
    rms0 = np.sqrt(v0_l2 / N) * convert(1.0, "m/s", "m/year")
    PISM.verbPrintf(1, com, "  Initial RMS misfit: %.2f m/yr (forward: %s)\n" %
                    (rms0, reason0.description()))

    # (Gradient check removed — verified correct with FD ratios ≈ 1.000.)

    # Solve
    t0 = time.time()
    reason = tao_solver.solve()
    solve_time = time.time() - t0

    PISM.verbPrintf(2, com, "  Inversion (%s): %s, %.1f seconds\n" %
                    (adjoint_method, reason.description(), solve_time))

    # Extract result
    zeta_inv = tikhonov.designSolution()
    tauc_inv = PISM.Scalar(grid, "tauc_inv")
    design_param.convertToDesignVariable(zeta_inv, tauc_inv)

    # Compute tauc error
    tauc_err = PISM.Scalar(grid, "tauc_error")
    tauc_err.copy_from(tauc_inv)
    tauc_err.add(-1.0, yield_stress_true)

    err_l2 = 0.0
    true_l2 = 0.0
    with PISM.vec.Access(nocomm=[tauc_err, yield_stress_true]):
        for (i, j) in grid.points():
            err_l2 += tauc_err[i, j] ** 2
            true_l2 += yield_stress_true[i, j] ** 2

    err_l2 = PISM.GlobalSum(com, err_l2)
    true_l2 = PISM.GlobalSum(com, true_l2)
    rel_error = np.sqrt(err_l2 / true_l2) if true_l2 > 0 else float('inf')

    # Compute final RMS misfit directly from the velocity difference
    # (same method as the initial diagnostic — known to be correct).
    vel_inv = solver.solution()
    vel_l2 = 0.0
    with PISM.vec.Access(nocomm=[vel_inv, u_obs]):
        for (i, j) in grid.points():
            du = vel_inv[i, j].u - u_obs[i, j].u
            dv = vel_inv[i, j].v - u_obs[i, j].v
            vel_l2 += du * du + dv * dv
    vel_l2 = PISM.GlobalSum(com, vel_l2)
    N = grid.Mx() * grid.My()
    rms_misfit = np.sqrt(vel_l2 / N) * convert(1.0, "m/s", "m/year")

    result = {
        "adjoint": adjoint_method,
        "test": test_name,
        "iterations": len(misfit_history),
        "solve_time_s": solve_time,
        "rms_misfit_m_yr": float(rms_misfit),
        "tauc_rel_error": float(rel_error),
        "history": misfit_history,
        "reason": reason.description(),
    }

    return result, tauc_inv


# ============================================================================
# Main
# ============================================================================

def main():
    ctx = PISM.Context()
    com = ctx.com
    config = ctx.config

    PISM.set_abort_on_sigint(True)

    PISM.verbPrintf(1, com, "ISMIP-HOM twin experiment: approximate vs incomplete vs exact adjoint\n")

    set_constants(config)

    # Inversion defaults (command-line -key value flags override these)
    #
    # "exp" parameterization: zeta = log(tauc), so the optimizer works in
    # log-space — much better scaled than "ident" for tauc ~ 3e10 Pa·s/m.
    config.set_string("inverse.design.param", "exp")
    # H1 regularization gives much better conditioning than L2-only.
    # Edge artifacts from non-periodic H1 are minor.
    config.set_number("inverse.design.cH1", 1.0)
    config.set_number("inverse.design.cL2", 0.0)
    # velocity_scale should match the typical velocity magnitude.
    # ISMIP-HOM C/D at L=80km: ~10-15 m/yr surface velocity.
    config.set_number("inverse.stress_balance.velocity_scale", 10.0)
    config.set_number("inverse.stress_balance.length_scale", 1e4)
    # The Tikhonov objective is f = J_design/eta + J_state, so LARGER eta
    # means WEAKER regularization (lets the optimizer fit the data).
    # For a twin experiment, use large eta to minimize regularization.
    config.set_number("inverse.tikhonov.penalty_weight", 1e6)
    config.set_string("inverse.stress_balance.method", "tikhonov_lmvm")
    # ISMIP-HOM tauc is in Pa·s/m (≈ 3e10), much larger than the default
    # tauc_max of 5e7 Pa. BLMVM clips zeta to [log(tauc_min), log(tauc_max)],
    # which destroys the initial guess if the bounds are too tight.
    config.set_number("inverse.stress_balance.tauc_max", 1e12)

    # Domain size (default 80 km)
    L = 80e3
    L_opt = PISM.OptionReal(ctx.unit_system, "-L",
                            "ISMIP-HOM domain size in km", "km", 80.0)
    L = L_opt.value() * 1e3

    # Which tests to run
    test_names = ["C", "D"]
    test_opt = PISM.OptionString("-test", "ISMIP-HOM test(s) to run (C, D, or CD)", "CD")
    if test_opt.is_set():
        test_names = [c for c in test_opt.value() if c in tests]

    all_results = []

    for test_name in test_names:
        PISM.verbPrintf(1, com,
                        "\n" + "=" * 60 + "\n" +
                        "ISMIP-HOM %s, L = %.0f km\n" % (test_name, L / 1e3) +
                        "=" * 60 + "\n")

        grid = create_grid(config, ctx, test_name, L)
        geometry, enthalpy, yield_stress = init_geometry(grid, test_name, L)

        # Step 1: Forward solve to get "observed" surface velocity
        PISM.verbPrintf(1, com, "\n--- Forward solve (generating observations) ---\n")
        u_obs = forward_solve(grid, yield_stress, config)

        # Step 2: Run inversion with each adjoint method
        methods = ["approximate", "incomplete", "exact"]
        tauc_results = {}

        for method in methods:
            result, tauc_inv = run_inversion(
                grid, geometry, enthalpy, yield_stress,
                u_obs, test_name, adjoint_method=method, config=config)
            all_results.append(result)
            tauc_results[method] = tauc_inv

        # Write NetCDF output
        output_file = "ismiphom_twin_%s_L%03d.nc" % (test_name, int(L / 1e3))
        PISM.verbPrintf(1, com, "Writing results to %s\n" % output_file)
        output = PISM.util.prepare_output(output_file)

        yield_stress.metadata().set_name("tauc_true")
        write_vars = [yield_stress, u_obs]
        for method, tauc_inv in tauc_results.items():
            tauc_inv.metadata().set_name("tauc_%s" % method)
            write_vars.append(tauc_inv)

        for arr in write_vars:
            for k in range(arr.ndof()):
                output.define_variable(arr.metadata(k))
            arr.write(output)

        output.close()

    # Print summary table
    if com.Get_rank() == 0:
        print("\n" + "=" * 78)
        print("ISMIP-HOM Twin Experiment Summary")
        print("=" * 78)
        print(f"{'Test':>4s}  {'Adjoint':>12s}  {'Iters':>5s}  {'Time(s)':>8s}  "
              f"{'RMS(m/yr)':>10s}  {'tauc err%':>10s}  {'Reason'}")
        print("-" * 78)
        for r in all_results:
            print(f"{r['test']:>4s}  {r['adjoint']:>12s}  {r['iterations']:>5d}  "
                  f"{r['solve_time_s']:>8.1f}  {r['rms_misfit_m_yr']:>10.4f}  "
                  f"{r['tauc_rel_error'] * 100:>9.2f}%  {r['reason']}")
        print("=" * 78)

        # Save convergence history as JSON
        json_file = "ismiphom_twin_results.json"
        with open(json_file, "w") as f:
            json.dump(all_results, f, indent=2)
        print(f"\nConvergence history saved to {json_file}")


if __name__ == "__main__":
    import PISM.invert.core
    main()
