/* Copyright (C) 2024 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "pism/stressbalance/ssa/SSAFD_SNES.hh"
#include "pism/stressbalance/StressBalance.hh" // Inputs
#include "pism/util/petscwrappers/Vec.hh"
#include <algorithm>            // std::max()

namespace pism {
namespace stressbalance {

PetscErrorCode SSAFDSNESConvergenceTest(SNES snes, PetscInt it, PetscReal xnorm, PetscReal gnorm,
                                        PetscReal f, SNESConvergedReason *reason, void *ctx) {
  PetscErrorCode ierr;

  SSAFD_SNES *solver = reinterpret_cast<SSAFD_SNES *>(ctx);
  double tolerance = solver->tolerance();

  ierr = SNESConvergedDefault(snes, it, xnorm, gnorm, f, reason, ctx); CHKERRQ(ierr);
  if (*reason >= 0 and tolerance > 0) {
    // converged or iterating
    Vec residual;
    ierr = SNESGetFunction(snes, &residual, NULL, NULL);
    CHKERRQ(ierr);

    PetscReal norm;
    ierr = VecNorm(residual, NORM_INFINITY, &norm);
    CHKERRQ(ierr);

    if (norm <= tolerance) {
      *reason = SNES_CONVERGED_FNORM_ABS;
    }
  }

  return 0;
}

double SSAFD_SNES::tolerance() const {
  return m_config->get_number("stress_balance.ssa.fd.absolute_tolerance");
}

SSAFD_SNES::SSAFD_SNES(std::shared_ptr<const Grid> grid, bool regional_mode)
    : SSAFDBase(grid, regional_mode), m_residual(grid, "_ssa_residual") {

  PetscErrorCode ierr;

  int stencil_width=2;
  m_DA = m_grid->get_dm(2, stencil_width);

  // ierr = DMCreateGlobalVector(*m_DA, m_X.rawptr());
  // PISM_CHK(ierr, "DMCreateGlobalVector");

  ierr = SNESCreate(m_grid->com, m_snes.rawptr());
  PISM_CHK(ierr, "SNESCreate");

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  m_callback_data.da = *m_DA;
  m_callback_data.solver = this;
  m_callback_data.inputs = nullptr;

  ierr = DMDASNESSetFunctionLocal(*m_DA, INSERT_VALUES,
#if PETSC_VERSION_LT(3,21,0)
                                  (DMDASNESFunction)SSAFD_SNES::function_callback,
#else
                                  (DMDASNESFunctionFn*)SSAFD_SNES::function_callback,
#endif
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  ierr = DMDASNESSetJacobianLocal(*m_DA,
#if PETSC_VERSION_LT(3,21,0)
                                  (DMDASNESJacobian)SSAFD_SNES::jacobian_callback,
#else
                                  (DMDASNESJacobianFn*)SSAFD_SNES::jacobian_callback,
#endif
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetJacobianLocal");

  // ierr = DMSetMatType(*m_DA, "baij");
  // PISM_CHK(ierr, "DMSetMatType");

  ierr = DMSetApplicationContext(*m_DA, &m_callback_data);
  PISM_CHK(ierr, "DMSetApplicationContext");

  ierr = SNESSetOptionsPrefix(m_snes, "ssafd_");
  PISM_CHK(ierr, "SNESSetOptionsPrefix");

  ierr = SNESSetDM(m_snes, *m_DA);
  PISM_CHK(ierr, "SNESSetDM");

  ierr = SNESSetConvergenceTest(m_snes, SSAFDSNESConvergenceTest, this, NULL);
  PISM_CHK(ierr, "SNESSetConvergenceTest");

  ierr = SNESSetTolerances(m_snes, 0.0, 0.0, 0.0, 500, -1);
  PISM_CHK(ierr, "SNESSetTolerances");

  ierr = SNESSetFromOptions(m_snes);
  PISM_CHK(ierr, "SNESSetFromOptions");
}

void SSAFD_SNES::solve(const Inputs &inputs) {
  m_callback_data.inputs = &inputs;
  initialize_iterations(inputs);
  {
    PetscErrorCode ierr;

    // Solve:
    // ierr = SNESSolve(m_snes, NULL, m_X);
    ierr = SNESSolve(m_snes, NULL, m_velocity_global.vec());
    PISM_CHK(ierr, "SNESSolve");

    // See if it worked.
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(m_snes, &reason);
    PISM_CHK(ierr, "SNESGetConvergedReason");
    if (reason < 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "SSAFD_SNES solve failed to converge (SNES reason %s)",
                                    SNESConvergedReasons[reason]);
    }

    PetscInt snes_iterations = 0;
    ierr = SNESGetIterationNumber(m_snes, &snes_iterations);
    PISM_CHK(ierr, "SNESGetIterationNumber");

    PetscInt ksp_iterations = 0;
    ierr = SNESGetLinearSolveIterations(m_snes, &ksp_iterations);
    PISM_CHK(ierr, "SNESGetLinearSolveIterations");

    m_log->message(1, "SSA: %d*%d its, %s\n", (int)snes_iterations,
                   (int)(ksp_iterations / std::max((int)snes_iterations, 1)),
                   SNESConvergedReasons[reason]);
  }
  m_callback_data.inputs = nullptr;

  // copy from m_velocity_global to provide m_velocity with ghosts:
  m_velocity.copy_from(m_velocity_global);

  compute_residual(inputs, m_velocity, m_residual);
}


PetscErrorCode SSAFD_SNES::function_callback(DMDALocalInfo * /*unused*/,
                                             Vector2d const *const *velocity, Vector2d **result,
                                             CallbackData *data) {
  try {
    data->solver->compute_residual(*data->inputs, velocity, result);
  } catch (...) {
    MPI_Comm com        = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)data->da, &com);
    CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

void SSAFD_SNES::compute_jacobian(const Inputs &inputs, Vector2d const *const *const velocity,
                                  Mat J) {
  fd_operator(*inputs.geometry, inputs.bc_mask, m_bc_scaling, *inputs.basal_yield_stress,
              m_basal_sliding_law, velocity, m_nuH, m_cell_type, &J, nullptr);
}

PetscErrorCode SSAFD_SNES::jacobian_callback(DMDALocalInfo * /*unused*/,
                                             Vector2d const *const *const velocity, Mat /* A */,
                                             Mat J, CallbackData *data) {
  try {
    data->solver->compute_jacobian(*data->inputs, velocity, J);
  } catch (...) {
    MPI_Comm com        = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)data->da, &com);
    CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

const array::Vector &SSAFD_SNES::residual() const {
  return m_residual;
}

//! @brief Computes the magnitude of the driving shear stress at the base of
//! ice (diagnostically).
class SSAFD_residual_mag : public Diag<SSAFD_SNES> {
public:
  SSAFD_residual_mag(const SSAFD_SNES *m) : Diag<SSAFD_SNES>(m) {

    // set metadata:
    m_vars = { { m_sys, "ssa_residual_mag" } };

    m_vars[0].long_name("magnitude of the SSAFD solver's residual").units("Pa");
  }

protected:
  virtual std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Scalar>("ssa_residual_mag");
    result->metadata(0) = m_vars[0];

    compute_magnitude(model->residual(), *result);

    return result;
  }
};

DiagnosticList SSAFD_SNES::diagnostics_impl() const {
  DiagnosticList result = SSAFDBase::diagnostics_impl();

  result["ssa_residual"] = Diagnostic::wrap(m_residual);
  result["ssa_residual_mag"] = Diagnostic::Ptr(new SSAFD_residual_mag(this));

  return result;
}


} // namespace stressbalance
} // namespace pism
