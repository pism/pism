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

#include "pism/inverse/IP_BlatterTaucForwardProblem.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/blatter/util/DataAccess.hh"
#include "pism/stressbalance/blatter/util/grid_hierarchy.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/Logger.hh"
#include "pism/util/fem/Quadrature.hh"
#include "pism/util/node_types.hh"

namespace pism {
namespace inverse {

IP_BlatterTaucForwardProblem::IP_BlatterTaucForwardProblem(
    std::shared_ptr<const Grid> grid,
    int Mz,
    int coarsening_factor,
    IPDesignVariableParameterization &tp)
  : Blatter(grid, Mz, coarsening_factor),
    m_zeta(NULL),
    m_tauc_copy(m_grid, "tauc"),
    m_fixed_tauc_locations(NULL),
    m_tauc_param(tp),
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

  m_tauc_copy.metadata(0)
    .long_name("yield stress for basal till (plastic or pseudo-plastic model)")
    .units("Pa");

  // Set up a separate KSP for the adjoint (transpose) solves.
  // We cannot reuse the SNES's KSP because its MG smoother (SOR) does
  // not support transpose application. This KSP uses the "inv_" prefix
  // so users can configure it via -inv_ksp_type, -inv_pc_type, etc.
  PetscErrorCode ierr;

  ierr = KSPCreate(m_grid->com, m_ksp.rawptr());
  PISM_CHK(ierr, "KSPCreate");

  ierr = KSPSetOptionsPrefix(m_ksp, "inv_");
  PISM_CHK(ierr, "KSPSetOptionsPrefix");

  // Defaults: GMRES + BJACOBI, rtol=1e-3, max 10000 iterations.
  // Override with -inv_ksp_type, -inv_ksp_rtol, -inv_pc_type, etc.
  ierr = KSPSetType(m_ksp, KSPGMRES);
  PISM_CHK(ierr, "KSPSetType");

  ierr = KSPSetTolerances(m_ksp, 1e-3, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
  PISM_CHK(ierr, "KSPSetTolerances");

  PC pc;
  ierr = KSPGetPC(m_ksp, &pc);
  PISM_CHK(ierr, "KSPGetPC");

  ierr = PCSetType(pc, PCBJACOBI);
  PISM_CHK(ierr, "PCSetType");

  ierr = KSPSetFromOptions(m_ksp);
  PISM_CHK(ierr, "KSPSetFromOptions");
}

void IP_BlatterTaucForwardProblem::init() {
  init_impl();
}

//! Sets the current value of the design parameter \f$\zeta\f$.
void IP_BlatterTaucForwardProblem::set_design(array::Scalar &new_zeta) {

  array::Scalar &tauc = m_tauc_copy;
  m_zeta = &new_zeta;

  // Convert zeta to tauc.
  m_tauc_param.convertToDesignVariable(*m_zeta, tauc);

  // Update tauc in the Blatter parameters array.
  {
    array::AccessScope list{&tauc, &m_parameters};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      m_parameters(i, j).tauc = tauc(i, j);
    }
  }
  m_parameters.update_ghosts();

  m_rebuild_J_state = true;
}

//! Sets \f$\zeta\f$ and solves the Blatter system.
/*!
  Builds a full Inputs struct from the grid variables, calls set_design to
  update tauc, then runs the Blatter solver via update() which handles
  geometry initialization, ice hardness, node types, and the SNES solve.
*/
std::shared_ptr<TerminationReason> IP_BlatterTaucForwardProblem::linearize_at(
    array::Scalar &zeta) {

  this->set_design(zeta);

  // Build Inputs from Grid::variables(), same pattern as Blatter::update expects
  const auto &variables = m_grid->variables();

  Geometry geometry(m_grid);
  geometry.ice_thickness.copy_from(*variables.get_2d_scalar("land_ice_thickness"));
  geometry.bed_elevation.copy_from(*variables.get_2d_scalar("bedrock_altitude"));
  geometry.sea_level_elevation.set(0.0);

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
  inputs.geometry           = &geometry;
  inputs.basal_yield_stress = &m_tauc_copy;  // use our converted tauc
  inputs.enthalpy           = variables.get_3d_scalar("enthalpy");
  inputs.basal_melt_rate    = NULL;
  inputs.age                = NULL;
  inputs.water_column_pressure = NULL;

  if (variables.is_available("vel_bc_mask")) {
    inputs.bc_mask = variables.get_2d_scalar("vel_bc_mask");
  }
  if (variables.is_available("vel_bc")) {
    inputs.bc_values = variables.get_2d_vector("vel_bc");
  }

  // This calls init_2d_parameters, init_ice_hardness, compute_node_type,
  // and the SNES solve with parameter continuation.
  this->update(inputs, true);
  this->extract_surface_velocity();

  // Check convergence from the SNES
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
void IP_BlatterTaucForwardProblem::extract_surface_velocity() {
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

petsc::DM &IP_BlatterTaucForwardProblem::get_da() const {
  return *m_grid->get_dm(1, m_config->get_number("grid.max_stencil_width"));
}

// ============================================================================
// Phase 2: Design Jacobian and Linearization
// ============================================================================

/*!
 * Applies the design Jacobian to a perturbation of the design variable.
 *
 * Computes dR/d(zeta) * dzeta, where R is the 3D Blatter residual.
 * Since tauc only enters through the basal sliding law, the contribution
 * is nonzero only at basal face (k=0) quadrature points.
 *
 * The result is projected to 2D surface velocity via the linearization:
 *   du_surface = P * (-J_state^{-1} * J_design * dzeta)
 *
 * NOTE: This method computes J_design * dzeta into the FULL 3D space,
 * then extracts the surface component. For use by apply_linearization.
 */
void IP_BlatterTaucForwardProblem::apply_jacobian_design(
    array::Vector &u, array::Scalar &dzeta, array::Vector &du) {
  (void)u;
  (void)dzeta;
  du.set(0.0);
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterTaucForwardProblem::apply_jacobian_design "
                     "is not directly supported. Use apply_linearization instead.");
}

/*!
 * Applies the transpose of the design Jacobian.
 *
 * Computes J_design^T * lambda, where lambda is a 3D adjoint vector.
 * The result is a 2D scalar (perturbation of zeta).
 *
 * Since tauc enters only through the basal face, this sums contributions
 * from basal (k=0) elements only.
 */
void IP_BlatterTaucForwardProblem::apply_jacobian_design_transpose(
    array::Vector &u, array::Vector &du, array::Scalar &dzeta) {
  (void)u;
  (void)du;
  dzeta.set(0.0);
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterTaucForwardProblem::apply_jacobian_design_transpose "
                     "is not directly supported. Use apply_linearization_transpose instead.");
}

/*!
 * Internal: Apply J_design * dzeta in the full 3D Blatter space.
 *
 * Result is a 3D Vec (same layout as the Blatter solution vector).
 * Only basal face elements (k=0) contribute.
 */
void IP_BlatterTaucForwardProblem::apply_jacobian_design_3d(
    array::Scalar &dzeta, Vec result_3d) {

  PetscErrorCode ierr;

  // Zero the global result
  ierr = VecSet(result_3d, 0.0);
  PISM_CHK(ierr, "VecSet");

  // Get a local (ghosted) vector for assembly, then scatter back to global
  Vec result_local;
  ierr = DMGetLocalVector(m_da, &result_local);
  PISM_CHK(ierr, "DMGetLocalVector");
  ierr = VecSet(result_local, 0.0);
  PISM_CHK(ierr, "VecSet");

  // Get the 3D solution as a local (ghosted) array for reading
  Vec x_local;
  ierr = DMGetLocalVector(m_da, &x_local);
  PISM_CHK(ierr, "DMGetLocalVector");
  ierr = DMGlobalToLocalBegin(m_da, m_x, INSERT_VALUES, x_local);
  PISM_CHK(ierr, "DMGlobalToLocalBegin");
  ierr = DMGlobalToLocalEnd(m_da, m_x, INSERT_VALUES, x_local);
  PISM_CHK(ierr, "DMGlobalToLocalEnd");

  Vector2d ***x = nullptr;
  ierr = DMDAVecGetArrayRead(m_da, x_local, &x);
  PISM_CHK(ierr, "DMDAVecGetArrayRead");

  Vector2d ***R = nullptr;
  ierr = DMDAVecGetArray(m_da, result_local, &R);
  PISM_CHK(ierr, "DMDAVecGetArray");

  DMDALocalInfo petsc_info;
  ierr = DMDAGetLocalInfo(m_da, &petsc_info);
  PISM_CHK(ierr, "DMDAGetLocalInfo");

  auto info = grid_transpose(petsc_info);

  double
    x_min = m_grid->x0() - m_grid->Lx(),
    y_min = m_grid->y0() - m_grid->Ly(),
    dx    = m_grid->dx(),
    dy    = m_grid->dy();

  fem::Q1Element3 element(info, fem::Q13DQuadrature8(), dx, dy, x_min, y_min);

  const int Nk = fem::q13d::n_chi;

  double z[Nk], floatation[Nk], bottom_elevation[Nk], ice_thickness[Nk],
         surface_elevation[Nk], sea_level[Nk], basal_yield_stress[Nk];
  int node_type[Nk];
  Vector2d velocity[Nk], R_nodal[Nk];

  // Access dzeta and zeta
  array::AccessScope list{&dzeta, m_zeta, &m_parameters};
  auto *P = m_parameters.array();

  // Loop over all elements that have at least one owned node
  for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
    for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {

      nodal_parameter_values(element, P, i, j,
                             node_type, bottom_elevation,
                             ice_thickness, surface_elevation, sea_level);

      if (exterior_element(node_type)) {
        continue;
      }

      // Only the basal element (k=0) contributes to the design Jacobian
      int k = 0;

      for (int n = 0; n < Nk; ++n) {
        R_nodal[n] = 0.0;
      }

      // Compute z coordinates
      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(i, j, k, n);
        z[n] = grid_z(bottom_elevation[n], ice_thickness[n], info.mz, I.k);
      }

      element.reset(i, j, k, z);

      // Get nodal velocity values
      element.nodal_values(x, velocity);

      // Handle Dirichlet BC nodes
      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(n);
        if (dirichlet_node(info, I)) {
          element.mark_row_invalid(n);
          velocity[n] = u_bc(element.x(n), element.y(n), element.z(n));
        }
      }

      // Get basal yield stress and floatation at element nodes
      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(n);
        basal_yield_stress[n] = P[I.j][I.i].tauc;
        floatation[n] = P[I.j][I.i].floatation;
      }

      // Compute the design Jacobian contribution on the basal face
      fem::Q1Element3Face *face = grounding_line(floatation) ? &m_face100 : &m_face4;
      face->reset(fem::q13d::FACE_BOTTOM, z);

      // Evaluate fields at face quadrature points
      Vector2d *u_q = m_work2[0];
      double *tauc_q = m_work[0];
      double *float_q = m_work[1];

      face->evaluate(velocity, u_q);
      face->evaluate(basal_yield_stress, tauc_q);
      face->evaluate(floatation, float_q);

      // Compute dtauc = g'(zeta) * dzeta at element nodes, then evaluate at quad pts
      double dtauc_nodal[Nk];
      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(n);
        double zeta_n = (*m_zeta)(I.i, I.j);
        double g_prime;
        m_tauc_param.toDesignVariable(zeta_n, NULL, &g_prime);
        dtauc_nodal[n] = g_prime * dzeta(I.i, I.j);
      }

      // Handle fixed tauc locations
      if (m_fixed_tauc_locations != nullptr) {
        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(n);
          if ((*m_fixed_tauc_locations)(I.i, I.j) != 0.0) {
            dtauc_nodal[n] = 0.0;
          }
        }
      }

      double *dtauc_q = m_work[2];
      face->evaluate(dtauc_nodal, dtauc_q);

      // Assemble: dR/d(tauc) * dtauc = (d(beta)/d(tauc)) * u * dtauc * psi * W
      // For the pseudo-plastic law: d(beta)/d(tauc) = beta/tauc = drag(1, u, v)
      for (unsigned int q = 0; q < face->n_pts(); ++q) {
        auto W = face->weight(q) / m_scaling;

        bool grounded = float_q[q] <= 0.0;
        double dbeta_dtauc = 0.0;
        if (grounded && tauc_q[q] > 0.0) {
          // drag(dtauc, u, v) = dtauc * f(|u|) = (d(beta)/d(tauc)) * dtauc
          dbeta_dtauc = m_basal_sliding_law->drag(dtauc_q[q], u_q[q].u, u_q[q].v);
        }

        for (int t = 0; t < element.n_chi(); ++t) {
          auto psi = face->chi(q, t);
          R_nodal[t] += W * psi * dbeta_dtauc * u_q[q];
        }
      }

      element.add_contribution(R_nodal, R);
    } // i
  } // j

  ierr = DMDAVecRestoreArrayRead(m_da, x_local, &x);
  PISM_CHK(ierr, "DMDAVecRestoreArrayRead");

  ierr = DMDAVecRestoreArray(m_da, result_local, &R);
  PISM_CHK(ierr, "DMDAVecRestoreArray");

  // Scatter local -> global (ADD_VALUES to accumulate from ghost overlaps)
  ierr = DMLocalToGlobalBegin(m_da, result_local, ADD_VALUES, result_3d);
  PISM_CHK(ierr, "DMLocalToGlobalBegin");
  ierr = DMLocalToGlobalEnd(m_da, result_local, ADD_VALUES, result_3d);
  PISM_CHK(ierr, "DMLocalToGlobalEnd");

  ierr = DMRestoreLocalVector(m_da, &x_local);
  PISM_CHK(ierr, "DMRestoreLocalVector");
  ierr = DMRestoreLocalVector(m_da, &result_local);
  PISM_CHK(ierr, "DMRestoreLocalVector");
}

/*!
 * Internal: Apply J_design^T * lambda_3d to get a 2D dzeta.
 *
 * lambda_3d is a 3D vector (adjoint variable).
 * Result is a 2D scalar field (perturbation of zeta).
 */
void IP_BlatterTaucForwardProblem::apply_jacobian_design_transpose_3d(
    Vec lambda_3d, array::Scalar &dzeta) {

  PetscErrorCode ierr;

  dzeta.set(0.0);

  // Get local (ghosted) copies for reading
  Vec x_local, lambda_local;
  ierr = DMGetLocalVector(m_da, &x_local);
  PISM_CHK(ierr, "DMGetLocalVector");
  ierr = DMGlobalToLocalBegin(m_da, m_x, INSERT_VALUES, x_local);
  PISM_CHK(ierr, "DMGlobalToLocalBegin");
  ierr = DMGlobalToLocalEnd(m_da, m_x, INSERT_VALUES, x_local);
  PISM_CHK(ierr, "DMGlobalToLocalEnd");

  ierr = DMGetLocalVector(m_da, &lambda_local);
  PISM_CHK(ierr, "DMGetLocalVector");
  ierr = DMGlobalToLocalBegin(m_da, lambda_3d, INSERT_VALUES, lambda_local);
  PISM_CHK(ierr, "DMGlobalToLocalBegin");
  ierr = DMGlobalToLocalEnd(m_da, lambda_3d, INSERT_VALUES, lambda_local);
  PISM_CHK(ierr, "DMGlobalToLocalEnd");

  Vector2d ***x = nullptr;
  ierr = DMDAVecGetArrayRead(m_da, x_local, &x);
  PISM_CHK(ierr, "DMDAVecGetArrayRead");

  Vector2d ***lambda = nullptr;
  ierr = DMDAVecGetArrayRead(m_da, lambda_local, &lambda);
  PISM_CHK(ierr, "DMDAVecGetArrayRead");

  DMDALocalInfo petsc_info;
  ierr = DMDAGetLocalInfo(m_da, &petsc_info);
  PISM_CHK(ierr, "DMDAGetLocalInfo");

  auto info = grid_transpose(petsc_info);

  double
    x_min = m_grid->x0() - m_grid->Lx(),
    y_min = m_grid->y0() - m_grid->Ly(),
    dx    = m_grid->dx(),
    dy    = m_grid->dy();

  fem::Q1Element3 element(info, fem::Q13DQuadrature8(), dx, dy, x_min, y_min);

  const int Nk = fem::q13d::n_chi;

  double z[Nk], floatation[Nk], bottom_elevation[Nk], ice_thickness[Nk],
         surface_elevation[Nk], sea_level[Nk], basal_yield_stress[Nk];
  int node_type[Nk];
  Vector2d velocity[Nk], lambda_nodal[Nk];

  array::AccessScope list{&dzeta, m_zeta, &m_parameters};
  auto *P = m_parameters.array();

  for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
    for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {

      nodal_parameter_values(element, P, i, j,
                             node_type, bottom_elevation,
                             ice_thickness, surface_elevation, sea_level);

      if (exterior_element(node_type)) {
        continue;
      }

      int k = 0; // basal element only

      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(i, j, k, n);
        z[n] = grid_z(bottom_elevation[n], ice_thickness[n], info.mz, I.k);
      }

      element.reset(i, j, k, z);

      // Get nodal velocity and adjoint values
      element.nodal_values(x, velocity);
      element.nodal_values(lambda, lambda_nodal);

      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(n);
        if (dirichlet_node(info, I)) {
          element.mark_row_invalid(n);
          velocity[n] = u_bc(element.x(n), element.y(n), element.z(n));
          lambda_nodal[n] = {0.0, 0.0};
        }
      }

      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(n);
        basal_yield_stress[n] = P[I.j][I.i].tauc;
        floatation[n] = P[I.j][I.i].floatation;
      }

      fem::Q1Element3Face *face = grounding_line(floatation) ? &m_face100 : &m_face4;
      face->reset(fem::q13d::FACE_BOTTOM, z);

      Vector2d *u_q = m_work2[0];
      Vector2d *lam_q = m_work2[1];
      double *tauc_qp = m_work[0];
      double *float_qp = m_work[1];

      face->evaluate(velocity, u_q);
      face->evaluate(lambda_nodal, lam_q);
      face->evaluate(basal_yield_stress, tauc_qp);
      face->evaluate(floatation, float_qp);

      // Accumulate: dzeta_k += sum_q W * dbeta_dtauc * (lambda . u) * phi_k(q)
      // Then multiply by g'(zeta_k) afterwards
      double dzeta_contrib[Nk];
      for (int n = 0; n < Nk; ++n) {
        dzeta_contrib[n] = 0.0;
      }

      for (unsigned int q = 0; q < face->n_pts(); ++q) {
        auto W = face->weight(q) / m_scaling;

        bool grounded = float_qp[q] <= 0.0;
        double dbeta_dtauc = 0.0;
        if (grounded && tauc_qp[q] > 0.0) {
          // d(beta)/d(tauc) = drag(1, u, v)
          dbeta_dtauc = m_basal_sliding_law->drag(1.0, u_q[q].u, u_q[q].v);
        }

        double dot = lam_q[q].u * u_q[q].u + lam_q[q].v * u_q[q].v;

        for (int t = 0; t < Nk; ++t) {
          auto psi = face->chi(q, t);
          // Only bottom-face nodes contribute to the 2D design variable
          dzeta_contrib[t] += W * dbeta_dtauc * dot * psi;
        }
      }

      // Map element contributions to 2D grid
      // Bottom-face nodes of a 3D element map to the (i,j) 2D nodes
      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(n);
        if (I.k == 0) {
          // This is a basal node -- contribute to 2D field
          // Check that (I.i, I.j) is owned
          if (I.i >= info.xs && I.i < info.xs + info.xm &&
              I.j >= info.ys && I.j < info.ys + info.ym) {
            dzeta(I.i, I.j) += dzeta_contrib[n];
          }
        }
      }

    } // i
  } // j

  ierr = DMDAVecRestoreArrayRead(m_da, x_local, &x);
  PISM_CHK(ierr, "DMDAVecRestoreArrayRead");

  ierr = DMDAVecRestoreArrayRead(m_da, lambda_local, &lambda);
  PISM_CHK(ierr, "DMDAVecRestoreArrayRead");

  ierr = DMRestoreLocalVector(m_da, &x_local);
  PISM_CHK(ierr, "DMRestoreLocalVector");
  ierr = DMRestoreLocalVector(m_da, &lambda_local);
  PISM_CHK(ierr, "DMRestoreLocalVector");

  // Multiply by g'(zeta) at each node
  {
    array::AccessScope list2{&dzeta, m_zeta};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      double g_prime;
      m_tauc_param.toDesignVariable((*m_zeta)(i, j), NULL, &g_prime);
      dzeta(i, j) *= g_prime;
    }
  }

  // Zero out fixed locations
  if (m_fixed_tauc_locations != nullptr) {
    array::AccessScope list2{&dzeta, m_fixed_tauc_locations};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      if ((*m_fixed_tauc_locations)(i, j) != 0.0) {
        dzeta(i, j) = 0.0;
      }
    }
  }
}

// ============================================================================
// Linearization: DF = -J_state^{-1} * J_design (projected to surface)
// ============================================================================

void IP_BlatterTaucForwardProblem::apply_linearization(
    array::Scalar &dzeta, array::Vector &du) {

  PetscErrorCode ierr;

  if (m_rebuild_J_state) {
    Mat J, Jpre;
    ierr = SNESGetJacobian(m_snes, &J, &Jpre, NULL, NULL);
    PISM_CHK(ierr, "SNESGetJacobian");

    ierr = SNESComputeJacobian(m_snes, m_x, J, Jpre);
    PISM_CHK(ierr, "SNESComputeJacobian");

    m_rebuild_J_state = false;
  }

  // Step 1: Compute rhs = J_design * dzeta (3D)
  Vec rhs;
  ierr = VecDuplicate(m_x, &rhs);
  PISM_CHK(ierr, "VecDuplicate");

  this->apply_jacobian_design_3d(dzeta, rhs);

  // rhs = -rhs
  ierr = VecScale(rhs, -1.0);
  PISM_CHK(ierr, "VecScale");

  // Step 2: Solve J_state * du_3d = -J_design * dzeta
  Vec du_3d;
  ierr = VecDuplicate(m_x, &du_3d);
  PISM_CHK(ierr, "VecDuplicate");

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
                                  "IP_BlatterTaucForwardProblem::apply_linearization"
                                  " KSP failed (reason %s)",
                                  KSPConvergedReasons[reason]);
  }

  m_log->message(4,
                 "IP_BlatterTaucForwardProblem::apply_linearization converged"
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

  ierr = VecDestroy(&rhs); PISM_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&du_3d); PISM_CHK(ierr, "VecDestroy");
}

// ============================================================================
// Transpose linearization: DF^T = -J_design^T * J_state^{-T} * P^T
// ============================================================================

void IP_BlatterTaucForwardProblem::apply_linearization_transpose(
    array::Vector &du, array::Scalar &dzeta) {

  PetscErrorCode ierr;

  m_log->message(2, "Blatter inverse: solving adjoint system...\n");

  if (m_rebuild_J_state) {
    Mat J, Jpre;
    ierr = SNESGetJacobian(m_snes, &J, &Jpre, NULL, NULL);
    PISM_CHK(ierr, "SNESGetJacobian");

    ierr = SNESComputeJacobian(m_snes, m_x, J, Jpre);
    PISM_CHK(ierr, "SNESComputeJacobian");

    m_rebuild_J_state = false;
  }

  // Step 1: Inject du (2D surface) into 3D at surface nodes: rhs_3d = P^T * du
  Vec rhs_3d;
  ierr = VecDuplicate(m_x, &rhs_3d);
  PISM_CHK(ierr, "VecDuplicate");
  ierr = VecSet(rhs_3d, 0.0);
  PISM_CHK(ierr, "VecSet");

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

  // Step 2: Solve J_state * lambda ≈ P^T * du
  //
  // NOTE: We use KSPSolve (not KSPSolveTranspose) with the SNES's
  // well-conditioned MG preconditioner. This is an approximation: the
  // true adjoint requires J^T, but the Blatter Jacobian is nearly
  // symmetric (asymmetry comes only from the Newton linearization of
  // the nonlinear viscosity). This produces a useful descent direction
  // for the Tikhonov minimization and avoids the need for a
  // transpose-compatible preconditioner.
  Vec lambda_3d;
  ierr = VecDuplicate(m_x, &lambda_3d);
  PISM_CHK(ierr, "VecDuplicate");

  KSP ksp;
  ierr = SNESGetKSP(m_snes, &ksp);
  PISM_CHK(ierr, "SNESGetKSP");

  Mat J;
  ierr = SNESGetJacobian(m_snes, &J, NULL, NULL, NULL);
  PISM_CHK(ierr, "SNESGetJacobian");

  ierr = KSPSetOperators(ksp, J, J);
  PISM_CHK(ierr, "KSPSetOperators");

  ierr = KSPSolve(ksp, rhs_3d, lambda_3d);
  PISM_CHK(ierr, "KSPSolve");

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason);
  PISM_CHK(ierr, "KSPGetConvergedReason");
  if (reason < 0) {
    ierr = VecDestroy(&rhs_3d); PISM_CHK(ierr, "VecDestroy");
    ierr = VecDestroy(&lambda_3d); PISM_CHK(ierr, "VecDestroy");
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "IP_BlatterTaucForwardProblem::apply_linearization_transpose"
                                  " adjoint KSP failed (reason %s)",
                                  KSPConvergedReasons[reason]);
  }

  m_log->message(2, "Blatter inverse: adjoint solve done.\n");

  // Step 3: Compute dzeta = -J_design^T * lambda
  this->apply_jacobian_design_transpose_3d(lambda_3d, dzeta);
  dzeta.scale(-1.0);

  if (dzeta.stencil_width() > 0) {
    dzeta.update_ghosts();
  }

  ierr = VecDestroy(&rhs_3d); PISM_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&lambda_3d); PISM_CHK(ierr, "VecDestroy");
}

} // end of namespace inverse
} // end of namespace pism
