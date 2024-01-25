// Copyright (C) 2004--2024 Constantine Khroulev, Ed Bueler and Jed Brown
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cassert>
#include <stdexcept>

#include "pism/geometry/Geometry.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace stressbalance {

SSAFD::KSPFailure::KSPFailure(const char* reason)
  : RuntimeError(ErrorLocation(), std::string("SSAFD KSP (linear solver) failed: ") + reason){
  // empty
}

SSAFD::PicardFailure::PicardFailure(const std::string &message)
  : RuntimeError(ErrorLocation(), "SSAFD Picard iterations failed: " + message) {
  // empty
}

SSA* SSAFDFactory(std::shared_ptr<const Grid> g) {
  return new SSAFD(g);
}

/*!
Because the FD implementation of the SSA uses Picard iteration, a PETSc KSP
and Mat are used directly.  In particular we set up \f$A\f$
(Mat m_A) and a \f$b\f$ (= Vec m_rhs) and iteratively solve
linear systems
  \f[ A x = b \f]
where \f$x\f$ (= Vec m_velocity_global).  A PETSc SNES object is never created.
 */
SSAFD::SSAFD(std::shared_ptr<const Grid> grid)
    : SSA(grid),
      m_hardness(grid, "ice_hardness"),
      m_nuH(grid, "nuH"),
      m_nuH_old(grid, "nuH_old"),
      m_work(grid, "work_vector", array::WITH_GHOSTS, 1 /* stencil width */),
      m_cell_type(m_grid, "ssafd_cell_type"),
      m_rhs(grid, "right_hand_side"),
      m_taud(m_grid, "taud"),
      m_velocity_old(grid, "velocity_old"),
      m_scaling(1e9), // comparable to typical beta for an ice stream;
      m_Ax(grid, "matrix_times_solution") {

  m_velocity_old.metadata(0)
      .long_name("old SSA velocity field; used for re-trying with a different epsilon")
      .units("m s-1");

  m_hardness.metadata(0)
      .long_name("vertically-averaged ice hardness")
      .set_units_without_validation(pism::printf("Pa s^(1/%f)", m_flow_law->exponent()));

  m_nuH.metadata(0)
      .long_name("ice thickness times effective viscosity")
      .units("Pa s m");

  m_nuH_old.metadata(0)
      .long_name("ice thickness times effective viscosity (before an update)")
      .units("Pa s m");

  m_taud.metadata(0)
      .long_name("X-component of the driving shear stress at the base of ice")
      .units("Pa");
  m_taud.metadata(1)
      .long_name("Y-component of the driving shear stress at the base of ice")
      .units("Pa");

  m_work.metadata(0).long_name("temporary storage used to compute nuH");

  // grounded_dragging_floating integer mask
  m_cell_type.metadata(0)
      .long_name("ice-type (ice-free/grounded/floating/ocean) integer mask");
  m_cell_type.metadata()["flag_values"]   = { MASK_ICE_FREE_BEDROCK, MASK_GROUNDED, MASK_FLOATING,
                                         MASK_ICE_FREE_OCEAN };
  m_cell_type.metadata()["flag_meanings"] = "ice_free_bedrock grounded_ice floating_ice ice_free_ocean";

  // The nuH viewer:
  m_view_nuh        = false;
  m_nuh_viewer_size = 300;

  // PETSc objects and settings
  {
    auto dm = m_velocity_global.dm();

    PetscErrorCode ierr;
    ierr = DMSetMatType(*dm, MATAIJ);
    PISM_CHK(ierr, "DMSetMatType");

    ierr = DMCreateMatrix(*dm, m_A.rawptr());
    PISM_CHK(ierr, "DMCreateMatrix");

    ierr = KSPCreate(m_grid->com, m_KSP.rawptr());
    PISM_CHK(ierr, "KSPCreate");

    ierr = KSPSetOptionsPrefix(m_KSP, "ssafd_");
    PISM_CHK(ierr, "KSPSetOptionsPrefix");

    // Use non-zero initial guess (i.e. SSA velocities from the last
    // solve() call).
    ierr = KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);
    PISM_CHK(ierr, "KSPSetInitialGuessNonzero");

    // Use the initial residual norm.
    ierr = KSPConvergedDefaultSetUIRNorm(m_KSP);
    PISM_CHK(ierr, "KSPConvergedDefaultSetUIRNorm");
  }
}

//! @note Uses `PetscErrorCode` *intentionally*.
void SSAFD::pc_setup_bjacobi() {
  PetscErrorCode ierr;
  PC pc;

  ierr = KSPSetType(m_KSP, KSPGMRES);
  PISM_CHK(ierr, "KSPSetType");

  ierr = KSPSetOperators(m_KSP, m_A, m_A);
  PISM_CHK(ierr, "KSPSetOperators");

  // Get the PC from the KSP solver:
  ierr = KSPGetPC(m_KSP, &pc);
  PISM_CHK(ierr, "KSPGetPC");

  // Set the PC type:
  ierr = PCSetType(pc, PCBJACOBI);
  PISM_CHK(ierr, "PCSetType");

  // Process options:
  ierr = KSPSetFromOptions(m_KSP);
  PISM_CHK(ierr, "KSPSetFromOptions");
}

//! @note Uses `PetscErrorCode` *intentionally*.
void SSAFD::pc_setup_asm() {
  PetscErrorCode ierr;
  PC pc, sub_pc;

  // Set parameters equivalent to
  // -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu

  ierr = KSPSetType(m_KSP, KSPGMRES);
  PISM_CHK(ierr, "KSPSetType");

  ierr = KSPSetOperators(m_KSP, m_A, m_A);
  PISM_CHK(ierr, "KSPSetOperators");

  // Switch to using the "unpreconditioned" norm.
  ierr = KSPSetNormType(m_KSP, KSP_NORM_UNPRECONDITIONED);
  PISM_CHK(ierr, "KSPSetNormType");

  // Switch to "right" preconditioning.
  ierr = KSPSetPCSide(m_KSP, PC_RIGHT);
  PISM_CHK(ierr, "KSPSetPCSide");

  // Get the PC from the KSP solver:
  ierr = KSPGetPC(m_KSP, &pc);
  PISM_CHK(ierr, "KSPGetPC");

  // Set the PC type:
  ierr = PCSetType(pc, PCASM);
  PISM_CHK(ierr, "PCSetType");

  // Set the sub-KSP object to "preonly"
  KSP *sub_ksp;
  ierr = PCSetUp(pc);
  PISM_CHK(ierr, "PCSetUp");

  ierr = PCASMGetSubKSP(pc, NULL, NULL, &sub_ksp);
  PISM_CHK(ierr, "PCASMGetSubKSP");

  ierr = KSPSetType(*sub_ksp, KSPPREONLY);
  PISM_CHK(ierr, "KSPSetType");

  // Set the PC of the sub-KSP to "LU".
  ierr = KSPGetPC(*sub_ksp, &sub_pc);
  PISM_CHK(ierr, "KSPGetPC");

  ierr = PCSetType(sub_pc, PCLU);
  PISM_CHK(ierr, "PCSetType");

  // Let the user override all this:
  // Process options:
  ierr = KSPSetFromOptions(m_KSP);
  PISM_CHK(ierr, "KSPSetFromOptions");
}

void SSAFD::init_impl() {
  SSA::init_impl();

  m_log->message(2, "  [using the KSP-based finite difference implementation]\n");

  // options
  options::Integer viewer_size("-ssa_nuh_viewer_size", "nuH viewer size", m_nuh_viewer_size);
  m_nuh_viewer_size = viewer_size;
  m_view_nuh        = options::Bool("-ssa_view_nuh", "Enable the SSAFD nuH runtime viewer");

  if (m_config->get_flag("stress_balance.calving_front_stress_bc")) {
    m_log->message(2, "  using PISM-PIK calving-front stress boundary condition ...\n");
  }

  m_default_pc_failure_count     = 0;
  m_default_pc_failure_max_count = 5;
}

void SSAFD::assemble_matrix(const Inputs &inputs, const array::Vector1 &velocity,
                            const array::Staggered1 &nuH, const array::CellType1 &cell_type, Mat A) {
  array::AccessScope list{ &velocity };
  fd_operator(*inputs.geometry, inputs.bc_mask, *inputs.basal_yield_stress, velocity.array(), nuH,
              cell_type, &A, nullptr);
}


//! \brief Compute the vertically-averaged horizontal velocity from the shallow
//! shelf approximation.
/*!
This is the main procedure in the SSAFD.  It manages the nonlinear solve process
and the Picard iteration.

The outer loop (over index `k`) is the nonlinear iteration.  In this loop the effective
viscosity is computed by compute_nuH() and then the linear system is
set up and solved.

Specifically, we call the PETSc procedure KSPSolve() to solve the linear system.
Solving the linear system is also a loop, an iteration, but it occurs
inside KSPSolve().  The user has full control of the KSP solve through the PETSc
interface.  The default choicess for KSP type `-ksp_type` and preconditioner type
`-pc_type` are GMRES(30) for the former and block Jacobi with ILU on the
blocks for the latter.  The defaults usually work.  These choices are important
but poorly understood.  The eigenvalues of the linearized
SSA vary with ice sheet geometry and temperature in ways that are not well-studied.
Nonetheless these eigenvalues determine the convergence of
this (inner) linear iteration.  A well-chosen preconditioner can put the eigenvalues
in the right place so that the KSP can converge quickly.

The preconditioner will behave differently on different numbers of
processors.  If the user wants the results of SSA calculations to be
independent of the number of processors, then `-pc_type none` could
be used, but performance will be poor.

If you want to test different KSP methods, it may be helpful to see how many
iterations were necessary.  Use `-ksp_monitor`.
Initial testing implies that CGS takes roughly half the iterations of
GMRES(30), but is not significantly faster because the iterations are each
roughly twice as slow.  The outputs of PETSc options `-ksp_monitor_singular_value`,
`-ksp_compute_eigenvalues` and `-ksp_plot_eigenvalues -draw_pause N`
(the last holds plots for N seconds) may be useful to diagnose.

The outer loop terminates when the effective viscosity times thickness
is no longer changing much, according to the tolerance set by the
option `-ssafd_picard_rtol`. The outer loop also terminates when a maximum
number of iterations is exceeded. We save the velocity from the last
time step in order to have a better estimate of the effective
viscosity than the u=v=0 result.

In truth there is an "outer outer" loop (over index `l`).  This attempts to
over-regularize the effective viscosity if the nonlinear iteration (the "outer"
loop over `k`) is not converging with the default regularization.  The same
over-regularization is attempted if the KSP object reports that it has not
converged.

(An alternative for recovery in the KSP diverged case, suggested by Jed, is to
revert to a direct linear solve, either for the whole domain (not scalable) or
on the subdomains.  This recovery alternative requires a more nontrivial choice
but it may be worthwhile.  Note the user can already do `-pc_type asm
-sub_pc_type lu` at the command line, forcing subdomain direct solves.)

FIXME: update this doxygen comment
*/
void SSAFD::solve(const Inputs &inputs) {

  // These computations do not depend on the solution, so they need to
  // be done only once.
  {
    initialize_iterations(inputs);
    assemble_rhs(inputs, m_cell_type, m_taud, m_rhs);
  }

  // Store away old SSA velocity (it might be needed in case a solver
  // fails).
  m_velocity_old.copy_from(m_velocity);

  for (unsigned int k = 0; k < 3; ++k) {
    try {
      if (k == 0) {
        // default strategy
        picard_iteration(inputs, m_config->get_number("stress_balance.ssa.epsilon"), 1.0);

        break;
      }
      if (k == 1) {
        // try underrelaxing the iteration
        const double underrelax =
            m_config->get_number("stress_balance.ssa.fd.nuH_iter_failure_underrelaxation");
        m_log->message(
            1, "  re-trying with effective viscosity under-relaxation (parameter = %.2f) ...\n",
            underrelax);
        picard_iteration(inputs, m_config->get_number("stress_balance.ssa.epsilon"), underrelax);

        break;
      }
      if (k == 2) {
        // try over-regularization
        picard_strategy_regularization(inputs);

        break;
      }

      // if we reached this, then all strategies above failed
      write_system_petsc("all_strategies_failed");
      throw RuntimeError(PISM_ERROR_LOCATION, "all SSAFD strategies failed");
    } catch (PicardFailure &f) {
      // proceed to the next strategy
    }
  }

  if (m_config->get_flag("stress_balance.ssa.fd.extrapolate_at_margins")) {
    extrapolate_velocity(inputs.geometry->cell_type, m_velocity);
  }

  // Post-process velocities if the user asked for it:
  if (m_config->get_flag("stress_balance.ssa.fd.brutal_sliding")) {
    const double brutal_sliding_scaleFactor =
        m_config->get_number("stress_balance.ssa.fd.brutal_sliding_scale");
    m_velocity.scale(brutal_sliding_scaleFactor);

    m_velocity.update_ghosts();
  }
}

void SSAFD::picard_iteration(const Inputs &inputs, double nuH_regularization,
                             double nuH_iter_failure_underrelax) {

  if (m_default_pc_failure_count < m_default_pc_failure_max_count) {
    // Give BJACOBI another shot if we haven't tried it enough yet

    try {
      pc_setup_bjacobi();
      picard_manager(inputs, nuH_regularization, nuH_iter_failure_underrelax);

    } catch (KSPFailure &f) {

      m_default_pc_failure_count += 1;

      m_log->message(1, "  re-trying using the Additive Schwarz preconditioner...\n");

      pc_setup_asm();

      m_velocity.copy_from(m_velocity_old);

      picard_manager(inputs, nuH_regularization, nuH_iter_failure_underrelax);
    }

  } else {
    // otherwise use ASM
    pc_setup_asm();

    picard_manager(inputs, nuH_regularization, nuH_iter_failure_underrelax);
  }
}

//! \brief Manages the Picard iteration loop.
void SSAFD::picard_manager(const Inputs &inputs, double nuH_regularization,
                           double nuH_iter_failure_underrelax) {
  PetscErrorCode ierr;
  // ksp_iterations should be a PetscInt because it is used in the
  // KSPGetIterationNumber() call below
  PetscInt ksp_iterations, ksp_iterations_total = 0, outer_iterations;
  KSPConvergedReason reason;

  int max_iterations =
      static_cast<int>(m_config->get_number("stress_balance.ssa.fd.max_iterations"));
  double ssa_relative_tolerance =
      m_config->get_number("stress_balance.ssa.fd.relative_convergence");
  bool verbose = m_log->get_threshold() >= 2, very_verbose = m_log->get_threshold() > 2;

  // set the initial guess:
  m_velocity_global.copy_from(m_velocity);

  m_stdout_ssa.clear();

  {
    array::AccessScope scope{ &m_velocity };
    compute_nuH(inputs.geometry->ice_thickness, m_cell_type, m_velocity.array(), m_hardness,
                nuH_regularization, m_nuH);
  }
  if (m_view_nuh) {
    update_nuH_viewers(m_nuH);
  }

  // outer loop
  for (int k = 0; k < max_iterations; ++k) {

    if (very_verbose) {
      m_stdout_ssa += pism::printf("  %2d:", k);
    }

    // assemble (or re-assemble) matrix, which depends on updated viscosity
    assemble_matrix(inputs, m_velocity, m_nuH, m_cell_type, m_A);

    if (very_verbose) {
      m_stdout_ssa += "A:";
    }

    // Call PETSc to solve linear system by iterative method; "inner iteration":
    ierr = KSPSetOperators(m_KSP, m_A, m_A);
    PISM_CHK(ierr, "KSPSetOperator");

    ierr = KSPSolve(m_KSP, m_rhs.vec(), m_velocity_global.vec());
    PISM_CHK(ierr, "KSPSolve");

    // Check if diverged; report to standard out about iteration
    ierr = KSPGetConvergedReason(m_KSP, &reason);
    PISM_CHK(ierr, "KSPGetConvergedReason");

    if (reason < 0) {
      // KSP diverged
      m_log->message(1, "PISM WARNING:  KSPSolve() reports 'diverged'; reason = %d = '%s'\n",
                     reason, KSPConvergedReasons[reason]);

      write_system_petsc("kspdivergederror");

      // Tell the caller that we failed. (The caller might try again,
      // though.)
      throw KSPFailure(KSPConvergedReasons[reason]);
    }

    // report on KSP success; the "inner" iteration is done
    ierr = KSPGetIterationNumber(m_KSP, &ksp_iterations);
    PISM_CHK(ierr, "KSPGetIterationNumber");

    ksp_iterations_total += ksp_iterations;

    if (very_verbose) {
      m_stdout_ssa += pism::printf("S:%d,%d: ", (int)ksp_iterations, reason);
    }

    // limit ice speed
    {
      auto max_speed = m_config->get_number("stress_balance.ssa.fd.max_speed", "m second-1");
      int high_speed_counter = 0;

      array::AccessScope list{ &m_velocity_global };

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();

        auto speed = m_velocity_global(i, j).magnitude();

        if (speed > max_speed) {
          m_velocity_global(i, j) *= max_speed / speed;
          high_speed_counter += 1;
        }
      }

      high_speed_counter = GlobalSum(m_grid->com, high_speed_counter);

      if (high_speed_counter > 0) {
        m_log->message(2, "  SSA speed was capped at %d locations\n", high_speed_counter);
      }
    }

    // Communicate so that we have stencil width for evaluation of effective
    // viscosity on next "outer" iteration (and geometry etc. if done):
    // Note that copy_from() updates ghosts of m_velocity.
    m_velocity.copy_from(m_velocity_global);

    // update viscosity
    double nuH_norm = 0.0, nuH_norm_change = 0.0;
    {
      // in preparation of measuring change of effective viscosity:
      m_nuH_old.copy_from(m_nuH);

      {
        array::AccessScope scope{ &m_velocity };
        compute_nuH(inputs.geometry->ice_thickness, m_cell_type, m_velocity.array(), m_hardness,
                    nuH_regularization, m_nuH);
      }

      if (nuH_iter_failure_underrelax != 1.0) {
        m_nuH.scale(nuH_iter_failure_underrelax);
        m_nuH.add(1.0 - nuH_iter_failure_underrelax, m_nuH_old);
      }
      auto norm       = compute_nuH_norm(m_nuH, m_nuH_old);
      nuH_norm        = norm[0];
      nuH_norm_change = norm[1];
    }

    if (m_view_nuh) {
      update_nuH_viewers(m_nuH);
    }

    if (very_verbose) {
      m_stdout_ssa += pism::printf("|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n", nuH_norm,
                                   nuH_norm_change / nuH_norm);

      // assume that high verbosity shows interest in immediate
      // feedback about SSA iterations
      m_log->message(2, m_stdout_ssa);

      m_stdout_ssa.clear();
    }

    outer_iterations = k + 1;

    // check for viscosity convergence
    if (nuH_norm == 0 || nuH_norm_change / nuH_norm < ssa_relative_tolerance) {
      goto done;
    }
  } // outer loop (k)

  // If we're here, it means that we exceeded max_iterations and still
  // failed.

  throw PicardFailure(pism::printf("effective viscosity not converged after %d iterations\n"
                                   "with nuH_regularization=%8.2e.",
                                   max_iterations, nuH_regularization));

done:

  if (very_verbose) {
    auto tempstr =
        pism::printf("... =%5d outer iterations, ~%3.1f KSP iterations each\n",
                     (int)outer_iterations, ((double)ksp_iterations_total) / outer_iterations);
    m_stdout_ssa += tempstr;
  } else if (verbose) {
    // at default verbosity, just record last nuH_norm_change and iterations
    auto tempstr =
        pism::printf("%5d outer iterations, ~%3.1f KSP iterations each\n", (int)outer_iterations,
                     ((double)ksp_iterations_total) / outer_iterations);

    m_stdout_ssa += tempstr;
  }

  if (verbose) {
    m_stdout_ssa = "  SSA: " + m_stdout_ssa;
  }
}

//! Old SSAFD recovery strategy: increase the SSA regularization parameter.
void SSAFD::picard_strategy_regularization(const Inputs &inputs) {
  // this has no units; epsilon goes up by this ratio when previous value failed
  const double DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;
  double nuH_regularization                   = m_config->get_number("stress_balance.ssa.epsilon");
  unsigned int k = 0, max_tries = 5;

  if (nuH_regularization <= 0.0) {
    throw PicardFailure("will not attempt over-regularization of nuH\n"
                        "because nuH_regularization == 0.0.");
  }

  while (k < max_tries) {
    m_velocity.copy_from(m_velocity_old);
    m_log->message(1, "  re-trying with nuH_regularization multiplied by %8.2f...\n",
                   DEFAULT_EPSILON_MULTIPLIER_SSA);

    nuH_regularization *= DEFAULT_EPSILON_MULTIPLIER_SSA;

    try {
      // 1.0 is the under-relaxation parameter
      picard_iteration(inputs, nuH_regularization, 1.0);
      // if this call succeeded, stop over-regularizing
      break;
    } catch (PicardFailure &f) {
      k += 1;

      if (k == max_tries) {
        throw PicardFailure("exceeded the max. number of nuH over-regularization attempts");
      }
    }
  }
}

//! \brief Compute the norm of `nu H` and the norm of the change in `nu H`.
/*!
Verification and PST experiments
suggest that an \f$L^1\f$ criterion for convergence is best.  For verification
there seems to be little difference, presumably because the solutions are smooth
and the norms are roughly equivalent on a subspace of smooth functions.  For PST,
the \f$L^1\f$ criterion gives faster runs with essentially the same results.
Presumably that is because rapid (temporal and spatial) variation in
\f$\nu H\f$ occurs at margins, occupying very few horizontal grid cells.
For the significant (e.g.~in terms of flux) parts of the flow, it is o.k. to ignore
a bit of bad behavior at these few places, and \f$L^1\f$ ignores it more than
\f$L^2\f$ (much less \f$L^\infty\f$, which might not work at all).

Note: clobbers `nuH_old` -- it is used to compute the change in `nu H`.
 */
std::array<double, 2> SSAFD::compute_nuH_norm(const array::Staggered &nuH,
                                              array::Staggered &nuH_old) {

  const double area      = m_grid->cell_area();
  const NormType MY_NORM = NORM_1;

  // Test for change in nu
  nuH_old.add(-1, nuH);

  std::vector<double> nuNorm = nuH.norm(MY_NORM), nuChange = nuH_old.norm(MY_NORM);

  nuChange[0] *= area;
  nuChange[1] *= area;
  nuNorm[0] *= area;
  nuNorm[1] *= area;

  double norm_change = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));
  double norm        = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));

  return { norm, norm_change };
}


//! Update the nuH viewer, which shows log10(nu H).
void SSAFD::update_nuH_viewers(const array::Staggered &nuH) {

  array::Scalar tmp(m_grid, "nuH");
  tmp.metadata(0)
      .long_name("log10 of (viscosity * thickness)")
      .units("Pa s m");

  array::AccessScope list{&nuH, &tmp};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double avg_nuH = 0.5 * (nuH(i,j,0) + nuH(i,j,1));
    if (avg_nuH > 1.0e14) {
      tmp(i,j) = log10(avg_nuH);
    } else {
      tmp(i,j) = 14.0;
    }
  }

  if (not m_nuh_viewer) {
    m_nuh_viewer.reset(new petsc::Viewer(m_grid->com, "nuH", m_nuh_viewer_size,
                                         m_grid->Lx(), m_grid->Ly()));
  }

  tmp.view({m_nuh_viewer});
}

void SSAFD::write_system_petsc(const std::string &namepart) {
  PetscErrorCode ierr;

  // write a file with a fixed filename; avoid zillions of files
  std::string filename = "SSAFD_" + namepart + ".petsc";
  m_log->message(1,
             "  writing linear system to PETSc binary file %s ...\n", filename.c_str());

  petsc::Viewer viewer;       // will be destroyed automatically
  ierr = PetscViewerBinaryOpen(m_grid->com, filename.c_str(), FILE_MODE_WRITE,
                               viewer.rawptr());
  PISM_CHK(ierr, "PetscViewerBinaryOpen");

  ierr = MatView(m_A, viewer);
  PISM_CHK(ierr, "MatView");

  ierr = VecView(m_rhs.vec(), viewer);
  PISM_CHK(ierr, "VecView");
}

//! @brief Reports the nuH (viscosity times thickness) product on the staggered
//! grid.
class SSAFD_nuH : public Diag<SSAFD> {
public:
  SSAFD_nuH(const SSAFD *m) : Diag<SSAFD>(m) {
    m_vars = { { m_sys, "nuH[0]" }, { m_sys, "nuH[1]" } };
    m_vars[0]
        .long_name("ice thickness times effective viscosity, i-offset")
        .units("Pa s m")
        .output_units("kPa s m");
    m_vars[1]
        .long_name("ice thickness times effective viscosity, j-offset")
        .units("Pa s m")
        .output_units("kPa s m");
  }

protected:
  virtual std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Staggered>("nuH");

    result->copy_from(model->integrated_viscosity());

    return result;
  }
};

//! @brief Computes the driving shear stress at the base of ice
//! (diagnostically).
/*! This is *not* a duplicate of SSB_taud: SSAFD_taud::compute() uses
  SSAFD::compute_driving_stress(), which tries to be smarter near ice margins.
*/
class SSAFD_taud : public Diag<SSAFD> {
public:
  SSAFD_taud(const SSAFD *m) : Diag<SSAFD>(m) {

    // set metadata:
    m_vars = { { m_sys, "taud_x" }, { m_sys, "taud_y" } };

    m_vars[0].long_name("X-component of the driving shear stress at the base of ice");
    m_vars[1].long_name("Y-component of the driving shear stress at the base of ice");

    for (auto &v : m_vars) {
      v.units("Pa");
      v["comment"] = "this is the driving stress used by the SSAFD solver";
    }
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Vector>("taud");

    result->copy_from(model->driving_stress());

    return result;
  }
};


//! @brief Computes the magnitude of the driving shear stress at the base of
//! ice (diagnostically).
class SSAFD_taud_mag : public Diag<SSAFD> {
public:
  SSAFD_taud_mag(const SSAFD *m) : Diag<SSAFD>(m) {

    // set metadata:
    m_vars = { { m_sys, "taud_mag" } };

    m_vars[0].long_name("magnitude of the driving shear stress at the base of ice").units("Pa");
    m_vars[0]["comment"] = "this is the magnitude of the driving stress used by the SSAFD solver";
  }

protected:
  virtual std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Scalar>("taud_mag");

    compute_magnitude(model->driving_stress(), *result);

    return result;
  }
};

DiagnosticList SSAFD::diagnostics_impl() const {
  DiagnosticList result = ShallowStressBalance::diagnostics_impl();

  // replace these diagnostics
  result["taud"] = Diagnostic::Ptr(new SSAFD_taud(this));
  result["taud_mag"] = Diagnostic::Ptr(new SSAFD_taud_mag(this));
  result["nuH"] = Diagnostic::Ptr(new SSAFD_nuH(this));

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
