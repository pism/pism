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

#include "pism/inverse/IP_BlatterHardavForwardProblem.hh"
#include "pism/rheology/FlowLaw.hh"
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

IP_BlatterHardavForwardProblem::IP_BlatterHardavForwardProblem(
    std::shared_ptr<const Grid> grid,
    int Mz,
    int coarsening_factor,
    IPDesignVariableParameterization &tp)
  : IP_BlatterForwardProblem(grid, Mz, coarsening_factor, tp),
    m_hardav(m_grid, "hardav")
{
  // Same (non-validated) units string as the "hardav" diagnostic and
  // PISM.model.createAveragedHardnessVec(), so that reading/writing the
  // field does not trigger a units conversion.
  m_hardav.metadata(0)
    .long_name("vertically-averaged ice hardness (Glen B = A^{-1/n})")
    .set_units_without_validation("Pa s^(1/n)");
}

//! Sets the current value of the design parameter `zeta` and computes
//! hardav via the parameterization. The hardness is handed to the Blatter
//! solver through stressbalance::Inputs::averaged_hardness in linearize_at
//! (see fill_inputs).
void IP_BlatterHardavForwardProblem::set_design(array::Scalar &new_zeta) {

  this->store_zeta(new_zeta);

  // Convert zeta to hardav.
  m_design_param.convertToDesignVariable(*m_zeta, m_hardav);

  m_rebuild_J_state = true;
}

void IP_BlatterHardavForwardProblem::fill_inputs(stressbalance::Inputs &inputs) {
  if (inputs.basal_yield_stress == nullptr) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "IP_BlatterHardavForwardProblem requires the basal yield "
                       "stress 'tauc' in Grid::variables()");
  }

  // Blatter::init_ice_hardness replicates m_hardav(i,j) across the sigma
  // column instead of deriving hardness from enthalpy.
  inputs.averaged_hardness = &m_hardav;
}

// ============================================================================
// Design Jacobian — volume integral analog of IP_SSAHardavForwardProblem
// ============================================================================
//
// The viscous part of the Blatter residual (Blatter::residual_f) is
//
//   R_t = sum_q (W_q / m_scaling) * eta_q * F(u, psi_t),
//   eta_q = E * (1/2) B_q (eps + gamma_q)^p,   B_q = sum_n chi_n(q) B_n,
//
// which is linear in the nodal hardness B_n. With the column-constant design
// variable B(i,j) = g(zeta(i,j)) this gives
//
//   J_design * dzeta  =  residual_f(u, dB),   dB_n = g'(zeta) * dzeta at (I.i, I.j),
//
//   (J_design^T lambda)(i,j) = g'(zeta(i,j)) *
//       sum_{elements in column} sum_q (W_q / m_scaling) * (eta_q / B_q) *
//       [ l_x (4 u_x + 2 v_y) + l_y (u_y + v_x) + l_z u_z
//       + m_x (u_y + v_x) + m_y (2 u_x + 4 v_y) + m_z v_z ](q) * chi_n(q),
//
// summed over all 8 nodes n of every element in the column (this is the
// "integration over the 3D velocities": every sigma level contributes to the
// single 2D design value of its column). Here (l, m) are the two components
// of the 3D adjoint field and l_x = sum_t lambda_t dpsi_t/dx, etc.

/*!
 * Internal: Apply J_design * dzeta in the full 3D Blatter space.
 *
 * Result is a 3D Vec (same layout as the Blatter solution vector). All
 * elements of the ice column contribute. `dzeta` is ghosted and already
 * zeroed at fixed design locations (see IP_BlatterForwardProblem::apply_linearization).
 */
void IP_BlatterHardavForwardProblem::apply_jacobian_design_3d(
    array::Scalar &dzeta, Vec result_3d) {

  PetscErrorCode ierr;

  ierr = VecSet(result_3d, 0.0);
  PISM_CHK(ierr, "VecSet");

  // Local (ghosted) vector for assembly, scattered back to global at the end
  Vec result_local;
  ierr = DMGetLocalVector(m_da, &result_local);
  PISM_CHK(ierr, "DMGetLocalVector");
  ierr = VecSet(result_local, 0.0);
  PISM_CHK(ierr, "VecSet");

  // The 3D solution as a local (ghosted) array
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

  double z[Nk], bottom_elevation[Nk], ice_thickness[Nk],
         surface_elevation[Nk], sea_level[Nk], dB_nodal[Nk];
  int node_type[Nk];
  Vector2d velocity[Nk], R_nodal[Nk];

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

      // dB = g'(zeta) * dzeta is column-constant: same at every sigma level
      double dB_column[4];
      for (int n = 0; n < 4; ++n) {
        auto I = element.local_to_global(i, j, 0, n);
        double g_prime;
        m_design_param.toDesignVariable((*m_zeta)(I.i, I.j), NULL, &g_prime);
        dB_column[n] = g_prime * dzeta(I.i, I.j);
      }

      // loop over elements in the column
      for (int k = info.gzs; k < info.gzs + info.gzm - 1; k++) {

        for (int n = 0; n < Nk; ++n) {
          R_nodal[n] = 0.0;
        }

        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(i, j, k, n);
          z[n] = grid_z(bottom_elevation[n], ice_thickness[n], info.mz, I.k);
        }

        element.reset(i, j, k, z);

        element.nodal_values(x, velocity);

        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(n);
          if (dirichlet_node(info, I)) {
            element.mark_row_invalid(n);
            velocity[n] = u_bc(element.x(n), element.y(n), element.z(n));
          }
          // nodes n and n+4 share the same (i, j) column
          dB_nodal[n] = dB_column[n % 4];
        }

        // Linearity of eta in B: dR/dB * dB = residual_f evaluated with B = dB.
        // residual_f includes the enhancement factor and 1/m_scaling, exactly
        // as the state Jacobian does.
        residual_f(element, velocity, dB_nodal, R_nodal);

        element.add_contribution(R_nodal, R);
      } // k
    } // i
  } // j

  ierr = DMDAVecRestoreArrayRead(m_da, x_local, &x);
  PISM_CHK(ierr, "DMDAVecRestoreArrayRead");

  ierr = DMDAVecRestoreArray(m_da, result_local, &R);
  PISM_CHK(ierr, "DMDAVecRestoreArray");

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
 * lambda_3d is a 3D vector (adjoint variable). The result is a 2D scalar
 * field: every element in a column contributes to the design value of that
 * column.
 */
void IP_BlatterHardavForwardProblem::apply_jacobian_design_transpose_3d(
    Vec lambda_3d, array::Scalar &dzeta) {

  PetscErrorCode ierr;

  dzeta.set(0.0);

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

  double z[Nk], bottom_elevation[Nk], ice_thickness[Nk],
         surface_elevation[Nk], sea_level[Nk];
  int node_type[Nk];
  Vector2d velocity[Nk], lambda_nodal[Nk];

  // quadrature-point storage: velocity and its gradient, adjoint and its gradient
  Vector2d
    *u   = m_work2[0],
    *u_x = m_work2[1],
    *u_y = m_work2[2],
    *u_z = m_work2[3],
    *l   = m_work2[4],
    *l_x = m_work2[5],
    *l_y = m_work2[6],
    *l_z = m_work2[7];

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

      // Contributions to the four 2D nodes of this column, accumulated over
      // all elements (sigma levels) in the column.
      double column_contrib[4] = {0.0, 0.0, 0.0, 0.0};

      for (int k = info.gzs; k < info.gzs + info.gzm - 1; k++) {

        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(i, j, k, n);
          z[n] = grid_z(bottom_elevation[n], ice_thickness[n], info.mz, I.k);
        }

        element.reset(i, j, k, z);

        element.nodal_values(x, velocity);
        element.nodal_values(lambda, lambda_nodal);

        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(n);
          if (dirichlet_node(info, I)) {
            // Dirichlet rows of the residual do not depend on B
            velocity[n] = u_bc(element.x(n), element.y(n), element.z(n));
            lambda_nodal[n] = {0.0, 0.0};
          }
        }

        element.evaluate(velocity, u, u_x, u_y, u_z);
        element.evaluate(lambda_nodal, l, l_x, l_y, l_z);

        for (unsigned int q = 0; q < element.n_pts(); ++q) {
          auto W = element.weight(q) / m_scaling;

          double
            ux = u_x[q].u,
            uy = u_y[q].u,
            uz = u_z[q].u,
            vx = u_x[q].v,
            vy = u_y[q].v,
            vz = u_z[q].v;

          // second invariant, as in Blatter::residual_f
          double gamma = (ux * ux + vy * vy + ux * vy +
                          0.25 * ((uy + vx) * (uy + vx) + uz * uz + vz * vz));

          // eta / B (eta is linear in B), including the enhancement factor
          double eta_per_B;
          m_flow_law->effective_viscosity(1.0, gamma, m_viscosity_eps, &eta_per_B, nullptr);
          eta_per_B *= m_E_viscosity;

          // sum_t lambda_t . F(u, psi_t), using grad(lambda) = sum_t lambda_t grad(psi_t)
          double S = (l_x[q].u * (4.0 * ux + 2.0 * vy) +
                      l_y[q].u * (uy + vx) +
                      l_z[q].u * uz +
                      l_x[q].v * (uy + vx) +
                      l_y[q].v * (2.0 * ux + 4.0 * vy) +
                      l_z[q].v * vz);

          double c = W * eta_per_B * S;

          for (int n = 0; n < Nk; ++n) {
            // nodes n and n+4 share the same (i, j) column
            column_contrib[n % 4] += c * element.chi(q, n).val;
          }
        } // q
      } // k

      // Map column contributions to the 2D grid (owned nodes only)
      for (int n = 0; n < 4; ++n) {
        auto I = element.local_to_global(i, j, 0, n);
        if (I.i >= info.xs && I.i < info.xs + info.xm &&
            I.j >= info.ys && I.j < info.ys + info.ym) {
          dzeta(I.i, I.j) += column_contrib[n];
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
