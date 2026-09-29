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
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/blatter/util/DataAccess.hh"
#include "pism/stressbalance/blatter/util/grid_hierarchy.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/fem/Quadrature.hh"
#include "pism/util/node_types.hh"

namespace pism {
namespace inverse {

IP_BlatterTaucForwardProblem::IP_BlatterTaucForwardProblem(
    std::shared_ptr<const Grid> grid,
    int Mz,
    int coarsening_factor,
    IPDesignVariableParameterization &tp)
  : IP_BlatterForwardProblem(grid, Mz, coarsening_factor, tp),
    m_tauc_copy(m_grid, "tauc")
{
  m_tauc_copy.metadata(0)
    .long_name("yield stress for basal till (plastic or pseudo-plastic model)")
    .units("Pa");
}

//! Sets the current value of the design parameter \f$\zeta\f$.
void IP_BlatterTaucForwardProblem::set_design(array::Scalar &new_zeta) {

  this->store_zeta(new_zeta);

  array::Scalar &tauc = m_tauc_copy;

  // Convert zeta to tauc.
  m_design_param.convertToDesignVariable(*m_zeta, tauc);

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

//! Use our converted tauc; use an inverted `hardav` as ice hardness if one
//! is available (alternating tauc/hardav inversions).
void IP_BlatterTaucForwardProblem::fill_inputs(stressbalance::Inputs &inputs) {
  inputs.basal_yield_stress = &m_tauc_copy;

  const auto &variables = m_grid->variables();
  if (variables.is_available("hardav")) {
    inputs.averaged_hardness = variables.get_2d_scalar("hardav");
  }
}

// ============================================================================
// Design Jacobian (basal face integral)
// ============================================================================

/*!
 * Internal: Apply J_design * dzeta in the full 3D Blatter space.
 *
 * Result is a 3D Vec (same layout as the Blatter solution vector).
 * Only basal face elements (k=0) contribute. `dzeta` is ghosted and already
 * zeroed at fixed design locations (see IP_BlatterForwardProblem::apply_linearization).
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
        m_design_param.toDesignVariable(zeta_n, NULL, &g_prime);
        dtauc_nodal[n] = g_prime * dzeta(I.i, I.j);
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
      m_design_param.toDesignVariable((*m_zeta)(i, j), NULL, &g_prime);
      dzeta(i, j) *= g_prime;
    }
  }

  // Zero out fixed locations
  this->apply_fixed_locations(dzeta);
}

} // end of namespace inverse
} // end of namespace pism
