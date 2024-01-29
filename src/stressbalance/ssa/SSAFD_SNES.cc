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

namespace pism {
namespace stressbalance {

SSAFD_SNES::SSAFD_SNES(std::shared_ptr<const Grid> grid, bool regional_mode)
    : SSAFDBase(grid, regional_mode) {

  PetscErrorCode ierr;

  int stencil_width=2;
  m_DA = m_grid->get_dm(2, stencil_width);

  ierr = DMCreateGlobalVector(*m_DA, m_X.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");

  ierr = SNESCreate(m_grid->com, m_snes.rawptr());
  PISM_CHK(ierr, "SNESCreate");

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  m_callback_data.da = *m_DA;
  m_callback_data.solver = this;
  m_callback_data.inputs = nullptr;

  ierr = DMDASNESSetFunctionLocal(*m_DA, INSERT_VALUES,
                                  (DMDASNESFunction)SSAFD_SNES::function_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  ierr = DMDASNESSetJacobianLocal(*m_DA, (DMDASNESJacobian)SSAFD_SNES::jacobian_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetJacobianLocal");

  ierr = DMSetMatType(*m_DA, "baij");
  PISM_CHK(ierr, "DMSetMatType");

  ierr = DMSetApplicationContext(*m_DA, &m_callback_data);
  PISM_CHK(ierr, "DMSetApplicationContext");

  ierr = SNESSetOptionsPrefix(m_snes, "ssafd_snes");
  PISM_CHK(ierr, "SNESSetOptionsPrefix");

  ierr = SNESSetDM(m_snes, *m_DA);
  PISM_CHK(ierr, "SNESSetDM");

  ierr = SNESSetFromOptions(m_snes);
  PISM_CHK(ierr, "SNESSetFromOptions");
}

void SSAFD_SNES::solve(const Inputs &inputs) {
  m_callback_data.inputs = &inputs;
  initialize_iterations(inputs);
  {
    PetscErrorCode ierr;

    // Solve:
    ierr = SNESSolve(m_snes, NULL, m_X);
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

    m_log->message(1, "SSAFD_SNES converged (SNES reason %s)\n", SNESConvergedReasons[reason]);
  }
  m_callback_data.inputs = nullptr;

  // copy from m_velocity_global to provide m_velocity with ghosts:
  m_velocity.copy_from(m_velocity_global);
}


PetscErrorCode SSAFD_SNES::function_callback(DMDALocalInfo * /*unused*/, Vector2d const *const *velocity,
                                             Vector2d **f, CallbackData *data) {
  try {
    data->solver->compute_residual(*data->inputs, velocity, f);
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

PetscErrorCode SSAFD_SNES::jacobian_callback(DMDALocalInfo * /* info */,
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
} // namespace stressbalance
} // namespace pism
