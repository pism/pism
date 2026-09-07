// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#include "pism/inverse/IP_BlatterForwardProblem.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace inverse {

IP_BlatterForwardProblem::IP_BlatterForwardProblem(
    std::shared_ptr<const Grid> grid,
    int Mz,
    int coarsening_factor,
    IPDesignVariableParameterization &tp)
  : Blatter(grid, Mz, coarsening_factor),
    m_zeta(nullptr),
    m_zeta_local(m_grid, "zeta_local"),
    m_dzeta_local(m_grid, "dzeta_local"),
    m_fixed_design_locations(nullptr),
    m_design_param(tp),
    m_rebuild_J_state(true)
{
  m_surface_velocity.reset(new array::Vector(m_grid, "surface_velocity"));
  m_surface_velocity->metadata(0)
    .long_name("x-component of surface velocity from Blatter solver")
    .units("m s^-1")
    .output_units("m year^-1");
  m_surface_velocity->metadata(1)
    .long_name("y-component of surface velocity from Blatter solver")
    .units("m s^-1")
    .output_units("m year^-1");

  // Zero-initialize PETSc wrapper members. The Wrapper<T> class does NOT
  // zero-initialize m_value, so rawptr() always returns non-null (it's
  // &m_value), and get() returns garbage. We must explicitly null them.
  *m_J_picard.rawptr() = nullptr;
  *m_ksp.rawptr() = nullptr;
}

void IP_BlatterForwardProblem::init() {
  init_impl();
}

void IP_BlatterForwardProblem::store_zeta(array::Scalar &zeta) {
  if (zeta.stencil_width() > 0) {
    m_zeta = &zeta;
  } else {
    m_zeta_local.copy_from(zeta); // updates ghosts
    m_zeta = &m_zeta_local;
  }
}

void IP_BlatterForwardProblem::apply_fixed_locations(array::Scalar &dzeta) {
  if (m_fixed_design_locations == nullptr) {
    return;
  }

  array::AccessScope list{&dzeta, m_fixed_design_locations};
  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    if ((*m_fixed_design_locations)(i, j) != 0.0) {
      dzeta(i, j) = 0.0;
    }
  }
}

//! Sets \f$\zeta\f$ and solves the Blatter system.
/*!
  Builds a full Inputs struct from the grid variables (geometry, enthalpy,
  `tauc` when available, velocity boundary conditions), lets the derived
  class add the design variable via fill_inputs(), then runs the Blatter
  solver via update() which handles geometry initialization, ice hardness,
  node types, and the SNES solve.
*/
std::shared_ptr<TerminationReason> IP_BlatterForwardProblem::linearize_at(
    array::Scalar &zeta) {

  this->set_design(zeta);

  const auto &variables = m_grid->variables();

  Geometry geometry(m_grid);
  geometry.ice_thickness.copy_from(*variables.get_2d_scalar("land_ice_thickness"));
  geometry.bed_elevation.copy_from(*variables.get_2d_scalar("bedrock_altitude"));

  if (variables.is_available("sea_surface_height_above_reference_ellipsoid")) {
    geometry.sea_level_elevation.copy_from(
        *variables.get_2d_scalar("sea_surface_height_above_reference_ellipsoid"));
  } else {
    // Default: set sea level well below the bed to ensure all ice is grounded.
    geometry.sea_level_elevation.copy_from(geometry.bed_elevation);
    geometry.sea_level_elevation.shift(-1000.0);
  }

  if (m_config->get_flag("geometry.part_grid.enabled") &&
      variables.is_available("ice_area_specific_volume")) {
    geometry.ice_area_specific_volume.copy_from(
        *variables.get_2d_scalar("ice_area_specific_volume"));
  } else {
    geometry.ice_area_specific_volume.set(0.0);
  }

  geometry.ensure_consistency(
      m_config->get_number("stress_balance.ice_free_thickness_standard"));

  stressbalance::Inputs inputs;
  inputs.geometry              = &geometry;
  inputs.age                   = nullptr;
  inputs.water_column_pressure = nullptr;

  if (variables.is_available("tauc")) {
    inputs.basal_yield_stress = variables.get_2d_scalar("tauc");
  }
  if (variables.is_available("enthalpy")) {
    inputs.enthalpy = variables.get_3d_scalar("enthalpy");
  }
  if (variables.is_available("vel_bc_mask")) {
    inputs.bc_mask = variables.get_2d_scalar("vel_bc_mask");
  }
  if (variables.is_available("vel_bc")) {
    inputs.bc_values = variables.get_2d_vector("vel_bc");
  }

  // Design-variable-specific inputs (tauc copy, averaged hardness, ...).
  this->fill_inputs(inputs);

  if (inputs.basal_yield_stress == nullptr) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Blatter inverse forward problem: basal yield stress 'tauc' "
                       "is not available");
  }
  if (inputs.enthalpy == nullptr && inputs.averaged_hardness == nullptr) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Blatter inverse forward problem: neither 'enthalpy' nor "
                       "an averaged hardness is available to set ice hardness");
  }

  // This calls init_2d_parameters, init_ice_hardness, compute_node_type,
  // and the SNES solve with parameter continuation. A solver failure (e.g. a
  // line-search trial with an unphysical design variable) is reported as a
  // failed TerminationReason so that the minimizer can stop gracefully
  // instead of aborting the whole run.
  try {
    this->update(inputs, true);
  } catch (RuntimeError &e) {
    m_log->message(1, "Blatter inverse forward problem: forward solve failed: %s\n",
                   e.what());
    m_rebuild_J_state = true;
    return std::shared_ptr<TerminationReason>(
      new GenericTerminationReason(-1, "Blatter solve failed: " + std::string(e.what())));
  }
  this->extract_surface_velocity();

  // The state changed: the Newton Jacobian must be re-assembled at the new
  // solution before it is used for linearizations / adjoint solves.
  m_rebuild_J_state = true;

  SNESConvergedReason snes_reason;
  PetscErrorCode ierr = SNESGetConvergedReason(m_snes, &snes_reason);
  PISM_CHK(ierr, "SNESGetConvergedReason");

  if (snes_reason > 0) {
    return std::shared_ptr<TerminationReason>(
      new GenericTerminationReason(1, "Blatter solve converged"));
  } else {
    return std::shared_ptr<TerminationReason>(
      new GenericTerminationReason(-1, "Blatter solve failed to converge"));
  }
}

//! Extract the 2D surface velocity from the 3D Blatter solution.
void IP_BlatterForwardProblem::extract_surface_velocity() {
  PetscErrorCode ierr;
  Vector2d ***x = nullptr;

  ierr = DMDAVecGetArrayRead(m_da, m_x, &x);
  PISM_CHK(ierr, "DMDAVecGetArrayRead");

  int Mz = (int)m_u_sigma->levels().size();

  array::AccessScope list{m_surface_velocity.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    // Surface is at the top of the sigma grid
    (*m_surface_velocity)(i, j) = x[j][i][Mz - 1]; // STORAGE_ORDER
  }

  ierr = DMDAVecRestoreArrayRead(m_da, m_x, &x);
  PISM_CHK(ierr, "DMDAVecRestoreArrayRead");
}

petsc::DM &IP_BlatterForwardProblem::get_da() const {
  return *m_grid->get_dm(1, m_config->get_number("grid.max_stencil_width"));
}

void IP_BlatterForwardProblem::apply_jacobian_design(
    array::Vector &u, array::Scalar &dzeta, array::Vector &du) {
  (void)u;
  (void)dzeta;
  du.set(0.0);
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterForwardProblem::apply_jacobian_design "
                     "is not directly supported. Use apply_linearization instead.");
}

void IP_BlatterForwardProblem::apply_jacobian_design_transpose(
    array::Vector &u, array::Vector &du, array::Scalar &dzeta) {
  (void)u;
  (void)du;
  dzeta.set(0.0);
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterForwardProblem::apply_jacobian_design_transpose "
                     "is not directly supported. Use apply_linearization_transpose instead.");
}

//! Re-assemble the SNES Jacobian at the current solution m_x.
/*!
  After SNESSolve the matrix held by the SNES was assembled at the *previous*
  Newton iterate. Re-assembling at the converged solution makes the
  linearization and the adjoint consistent with the state that produced the
  surface velocity.
*/
void IP_BlatterForwardProblem::update_state_jacobian() {
  if (not m_rebuild_J_state) {
    return;
  }

  PetscErrorCode ierr;
  Mat J, Jpre;
  ierr = SNESGetJacobian(m_snes, &J, &Jpre, NULL, NULL);
  PISM_CHK(ierr, "SNESGetJacobian");

  ierr = SNESComputeJacobian(m_snes, m_x, J, Jpre);
  PISM_CHK(ierr, "SNESComputeJacobian");

  m_rebuild_J_state = false;
}

// ============================================================================
// Linearization: DF = -J_state^{-1} * J_design (projected to surface)
// ============================================================================

void IP_BlatterForwardProblem::apply_linearization(
    array::Scalar &dzeta, array::Vector &du) {

  PetscErrorCode ierr;

  this->update_state_jacobian();

  // Ghosted copy of dzeta with fixed design locations zeroed: the design
  // Jacobian kernels read it at ghost-element nodes.
  m_dzeta_local.copy_from(dzeta);
  this->apply_fixed_locations(m_dzeta_local);
  m_dzeta_local.update_ghosts();

  // Step 1: Compute rhs = J_design * dzeta (3D)
  Vec rhs;
  ierr = DMCreateGlobalVector(m_da, &rhs);
  PISM_CHK(ierr, "DMCreateGlobalVector");

  this->apply_jacobian_design_3d(m_dzeta_local, rhs);

  // rhs = -rhs
  ierr = VecScale(rhs, -1.0);
  PISM_CHK(ierr, "VecScale");

  // Step 2: Solve J_state * du_3d = -J_design * dzeta
  Vec du_3d;
  ierr = DMCreateGlobalVector(m_da, &du_3d);
  PISM_CHK(ierr, "DMCreateGlobalVector");

  // Reuse the SNES's KSP (configured via -bp_ksp_* flags)
  KSP ksp;
  ierr = SNESGetKSP(m_snes, &ksp);
  PISM_CHK(ierr, "SNESGetKSP");

  Mat J;
  ierr = SNESGetJacobian(m_snes, &J, NULL, NULL, NULL);
  PISM_CHK(ierr, "SNESGetJacobian");

  ierr = KSPSetOperators(ksp, J, J);
  PISM_CHK(ierr, "KSPSetOperators");

  ierr = KSPSolve(ksp, rhs, du_3d);
  PISM_CHK(ierr, "KSPSolve");

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason);
  PISM_CHK(ierr, "KSPGetConvergedReason");
  if (reason < 0) {
    ierr = VecDestroy(&rhs); PISM_CHK(ierr, "VecDestroy");
    ierr = VecDestroy(&du_3d); PISM_CHK(ierr, "VecDestroy");
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "IP_BlatterForwardProblem::apply_linearization"
                                  " KSP failed (reason %s)",
                                  KSPConvergedReasons[reason]);
  }

  m_log->message(4,
                 "IP_BlatterForwardProblem::apply_linearization converged"
                 " (KSP reason %s)\n",
                 KSPConvergedReasons[reason]);

  // Step 3: Extract surface velocity from du_3d
  {
    Vector2d ***du_arr = nullptr;
    ierr = DMDAVecGetArrayRead(m_da, du_3d, &du_arr);
    PISM_CHK(ierr, "DMDAVecGetArrayRead");

    int Mz = (int)m_u_sigma->levels().size();

    array::AccessScope list{&du};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      du(i, j) = du_arr[j][i][Mz - 1]; // surface
    }

    ierr = DMDAVecRestoreArrayRead(m_da, du_3d, &du_arr);
    PISM_CHK(ierr, "DMDAVecRestoreArrayRead");
  }

  if (du.stencil_width() > 0) {
    du.update_ghosts();
  }

  ierr = VecDestroy(&rhs); PISM_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&du_3d); PISM_CHK(ierr, "VecDestroy");
}

// ============================================================================
// Transpose linearization: DF^T = -J_design^T * J_state^{-T} * P^T
// ============================================================================

void IP_BlatterForwardProblem::apply_linearization_transpose(
    array::Vector &du, array::Scalar &dzeta) {

  PetscErrorCode ierr;

  this->update_state_jacobian();

  // Step 1: Create 3D vectors from the DM (not VecDuplicate) to ensure
  // proper DM association and parallel layout.
  Vec rhs_3d, lambda_3d;
  ierr = DMCreateGlobalVector(m_da, &rhs_3d);
  PISM_CHK(ierr, "DMCreateGlobalVector");
  ierr = DMCreateGlobalVector(m_da, &lambda_3d);
  PISM_CHK(ierr, "DMCreateGlobalVector");
  ierr = VecSet(rhs_3d, 0.0);
  PISM_CHK(ierr, "VecSet");

  // Inject du (2D surface) into 3D at surface nodes: rhs_3d = P^T * du
  {
    Vector2d ***rhs_arr = nullptr;
    ierr = DMDAVecGetArray(m_da, rhs_3d, &rhs_arr);
    PISM_CHK(ierr, "DMDAVecGetArray");

    int Mz = (int)m_u_sigma->levels().size();

    array::AccessScope list{&du};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      rhs_arr[j][i][Mz - 1] = du(i, j); // inject at surface
    }

    ierr = DMDAVecRestoreArray(m_da, rhs_3d, &rhs_arr);
    PISM_CHK(ierr, "DMDAVecRestoreArray");
  }

  // Standalone adjoint KSP (not the SNES's MG KSP).
  // Configure via -inv_adj_ksp_type, -inv_adj_pc_type, etc.
  if (m_ksp.get() == nullptr) {
    ierr = KSPCreate(m_grid->com, m_ksp.rawptr());
    PISM_CHK(ierr, "KSPCreate");

    ierr = KSPSetOptionsPrefix(m_ksp, "inv_adj_");
    PISM_CHK(ierr, "KSPSetOptionsPrefix");

    ierr = KSPSetType(m_ksp, KSPGMRES);
    PISM_CHK(ierr, "KSPSetType");

    ierr = KSPSetTolerances(m_ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
    PISM_CHK(ierr, "KSPSetTolerances");

    ierr = KSPSetFromOptions(m_ksp);
    PISM_CHK(ierr, "KSPSetFromOptions");
  }

  std::string adjoint_method = m_config->get_string("inverse.adjoint.method");

  if (adjoint_method == "incomplete") {
    // Incomplete (Picard) adjoint: assemble a separate Picard Jacobian
    // (drops viscosity derivative terms), then KSPSolve. The matrix is
    // truly symmetric, so CG + any preconditioner works.
    m_log->message(2, "Blatter inverse: adjoint solve (incomplete/Picard, KSPSolve)...\n");

    if (m_J_picard.get() == nullptr) {
      Mat J_snes;
      ierr = SNESGetJacobian(m_snes, &J_snes, NULL, NULL, NULL);
      PISM_CHK(ierr, "SNESGetJacobian");

      ierr = MatDuplicate(J_snes, MAT_DO_NOT_COPY_VALUES, m_J_picard.rawptr());
      PISM_CHK(ierr, "MatDuplicate");

      ierr = MatSetDM(m_J_picard, m_da);
      PISM_CHK(ierr, "MatSetDM");
    }

    this->compute_picard_jacobian(m_J_picard);

    ierr = KSPSetOperators(m_ksp, m_J_picard, m_J_picard);
    PISM_CHK(ierr, "KSPSetOperators");

    ierr = KSPSolve(m_ksp, rhs_3d, lambda_3d);
    PISM_CHK(ierr, "KSPSolve");

  } else if (adjoint_method == "approximate") {
    // Approximate adjoint: KSPSolve on the SNES Jacobian (symmetrized
    // by the upper-triangle mirror in compute_jacobian). Fast — reuses
    // the existing matrix, no reassembly. Any preconditioner works if
    // the matrix is approximately symmetric (GMRES recommended).
    m_log->message(2, "Blatter inverse: adjoint solve (approximate, KSPSolve)...\n");

    Mat J;
    ierr = SNESGetJacobian(m_snes, &J, NULL, NULL, NULL);
    PISM_CHK(ierr, "SNESGetJacobian");

    ierr = KSPSetOperators(m_ksp, J, J);
    PISM_CHK(ierr, "KSPSetOperators");

    ierr = KSPSolve(m_ksp, rhs_3d, lambda_3d);
    PISM_CHK(ierr, "KSPSolve");

  } else {
    // Exact adjoint: KSPSolveTranspose on the Newton Jacobian.
    // Requires transpose-compatible preconditioner (e.g., -inv_adj_pc_type jacobi).
    m_log->message(2, "Blatter inverse: adjoint solve (exact, KSPSolveTranspose)...\n");

    Mat J;
    ierr = SNESGetJacobian(m_snes, &J, NULL, NULL, NULL);
    PISM_CHK(ierr, "SNESGetJacobian");

    ierr = KSPSetOperators(m_ksp, J, J);
    PISM_CHK(ierr, "KSPSetOperators");

    ierr = KSPSolveTranspose(m_ksp, rhs_3d, lambda_3d);
    PISM_CHK(ierr, "KSPSolveTranspose");
  }

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(m_ksp, &reason);
  PISM_CHK(ierr, "KSPGetConvergedReason");

  PetscInt ksp_its;
  ierr = KSPGetIterationNumber(m_ksp, &ksp_its);
  PISM_CHK(ierr, "KSPGetIterationNumber");

  m_log->message(2, "  Adjoint KSP: %d iterations, reason: %s\n",
                 (int)ksp_its, KSPConvergedReasons[reason]);

  if (reason < 0) {
    ierr = VecDestroy(&rhs_3d); PISM_CHK(ierr, "VecDestroy");
    ierr = VecDestroy(&lambda_3d); PISM_CHK(ierr, "VecDestroy");
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "IP_BlatterForwardProblem::apply_linearization_transpose"
                                  " adjoint KSP failed (reason %s)",
                                  KSPConvergedReasons[reason]);
  }

  // Step 3: Compute dzeta = -J_design^T * lambda
  this->apply_jacobian_design_transpose_3d(lambda_3d, dzeta);
  dzeta.scale(-1.0);

  m_log->message(2, "Blatter inverse: adjoint solve done.\n");

  if (dzeta.stencil_width() > 0) {
    dzeta.update_ghosts();
  }

  ierr = VecDestroy(&rhs_3d); PISM_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&lambda_3d); PISM_CHK(ierr, "VecDestroy");
}

} // end of namespace inverse
} // end of namespace pism
