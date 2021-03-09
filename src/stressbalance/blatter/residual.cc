/* Copyright (C) 2020, 2021 PISM Authors
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

#include "Blatter.hh"

#include "pism/basalstrength/basal_resistance.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/util/node_types.hh"
#include "util/DataAccess.hh"
#include "util/grid_hierarchy.hh"    // grid_transpose(), grid_z()

namespace pism {
namespace stressbalance {

static const Vector2 u_exterior = {0.0, 0.0};

/*!
 * Computes the residual contribution of the "main" part of the Blatter system.
 */
void Blatter::residual_f(const fem::Q1Element3 &element,
                         const Vector2 *u_nodal,
                         const double *B_nodal,
                         Vector2 *residual) {

  Vector2
    *u   = m_work2[0],
    *u_x = m_work2[1],
    *u_y = m_work2[2],
    *u_z = m_work2[3];

  double *B = m_work[0];

  // evaluate u and its partial derivatives at quadrature points
  element.evaluate(u_nodal, u, u_x, u_y, u_z);

  // evaluate B (ice hardness) at quadrature points
  element.evaluate(B_nodal, B);

  // loop over all quadrature points
  for (int q = 0; q < element.n_pts(); ++q) {
    auto W = element.weight(q);

    double
      ux = u_x[q].u,
      uy = u_y[q].u,
      uz = u_z[q].u,
      vx = u_x[q].v,
      vy = u_y[q].v,
      vz = u_z[q].v;

    double gamma = (ux * ux + vy * vy + ux * vy +
                    0.25 * ((uy + vx) * (uy + vx) + uz * uz + vz * vz));

    double eta;
    m_flow_law->effective_viscosity(B[q], gamma, &eta, nullptr);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      const auto &psi = element.chi(q, t);

      residual[t].u += W * (eta * (psi.dx * (4.0 * ux + 2.0 * vy) +
                                   psi.dy * (uy + vx) +
                                   psi.dz * uz));
      residual[t].v += W * (eta * (psi.dx * (uy + vx) +
                                   psi.dy * (2.0 * ux + 4.0 * vy) +
                                   psi.dz * vz));
    }
  }
}

/*! Computes the residual contribution of the "source term".
 *
 * This term contains the driving stress.
 */
void Blatter::residual_source_term(const fem::Q1Element3 &element,
                                   const double *surface,
                                   const double *bed,
                                   Vector2 *residual) {

  // compute s_x and s_y
  double
    *s_x = m_work[0],
    *s_y = m_work[1];

  if (m_eta_transform) {
    // use the same storage for eta_x and eta_y
    double
      *eta_nodal = m_work[2],
      *eta       = m_work[3],
      *eta_x     = s_x,
      *eta_y     = s_y,
      *eta_z     = m_work[4];
    double
      *b   = m_work[5],
      *b_x = m_work[6],
      *b_y = m_work[7],
      *b_z = m_work[8];

    double
      n = m_glen_exponent,
      p = (2.0 * n + 2.0) / n;

    for (int k = 0; k < element.n_chi(); ++k) {
      eta_nodal[k] = pow(surface[k] - bed[k], p);
    }

    element.evaluate(eta_nodal, eta, eta_x, eta_y, eta_z);
    element.evaluate(bed, b, b_x, b_y, b_z);

    for (int q = 0; q < element.n_pts(); ++q) {
      double C = pow(eta[q], 1.0 / p - 1.0) / p;

      s_x[q] = C * eta_x[q] + b_x[q];
      s_y[q] = C * eta_y[q] + b_y[q];
    }
  } else {
    // these arrays are needed by the call below (but results are discarded)
    double
      *s   = m_work[2],
      *s_z = m_work[3];

    element.evaluate(surface, s, s_x, s_y, s_z);
  }

  for (int q = 0; q < element.n_pts(); ++q) {
    auto W = element.weight(q);

    auto F = m_rho_ice_g * Vector2(s_x[q], s_y[q]);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      const auto &psi = element.chi(q, t);

      residual[t] += W * psi.val * F;
    }
  }
}

/*!
 * Computes the residual contribution of the basal boundary condition.
 *
 * This takes care of basal sliding.
 */
void Blatter::residual_basal(const fem::Q1Element3 &element,
                             const fem::Q1Element3Face &face,
                             const double *tauc_nodal,
                             const double *f_nodal,
                             const Vector2 *u_nodal,
                             Vector2 *residual) {

  Vector2 *u = m_work2[0];

  double
    *tauc       = m_work[0],
    *floatation = m_work[1];

  face.evaluate(u_nodal, u);
  face.evaluate(tauc_nodal, tauc);
  face.evaluate(f_nodal, floatation);

  for (int q = 0; q < face.n_pts(); ++q) {
    auto W = face.weight(q);

    bool grounded = floatation[q] <= 0.0;
    double beta = grounded ? m_basal_sliding_law->drag(tauc[q], u[q].u, u[q].v) : 0.0;

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      residual[t].u += W * psi * beta * u[q].u;
      residual[t].v += W * psi * beta * u[q].v;
    }
  }
}

/*!
 * Top surface contribution to the residual.
 *
 * Used by verification tests ONLY.
 */
void Blatter::residual_surface(const fem::Q1Element3 &element,
                               const fem::Q1Element3Face &face,
                               Vector2 *residual) {
  (void) element;
  (void) face;
  (void) residual;
  // In normal circumstances the top surface contribution is zero (natural BCs apply).
}


/*!
 * Computes the residual contribution of lateral boundary conditions.
 *
 * This takes care of "calving front" stress boundary conditions.
 *
 * FIXME: make p_ocean an input from a parameterization of the melange back pressure.
 */
void Blatter::residual_lateral(const fem::Q1Element3 &element,
                               const fem::Q1Element3Face &face,
                               const double *surface_nodal,
                               const double *z_nodal,
                               const double *sl_nodal,
                               Vector2 *residual) {
  double
    *z         = m_work[0],
    *s         = m_work[1],
    *sea_level = m_work[2];

  // compute physical coordinates of quadrature points on this face
  face.evaluate(surface_nodal, s);
  face.evaluate(z_nodal, z);
  face.evaluate(sl_nodal, sea_level);

  // loop over all quadrature points
  for (int q = 0; q < face.n_pts(); ++q) {
    auto W = face.weight(q);
    auto N3 = face.normal(q);
    Vector2 N = {N3.x, N3.y};

    double
      ice_depth   = s[q] - z[q],
      p_ice       = m_rho_ice_g * ice_depth,
      water_depth = std::max(sea_level[q] - z[q], 0.0),
      p_ocean     = m_rho_ocean_g * water_depth;

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      residual[t] += - W * psi * (p_ice - p_ocean) * N;
    }
  }
}

/*! Set the residual at Dirichlet locations
 *
 * Compute the residual at Dirichlet locations and reset the residual to zero elsewhere.
 *
 * Setting it to zero is necessary because we call DMDASNESSetFunctionLocal() with
 * INSERT_VALUES.
 *
 */
void Blatter::residual_dirichlet(const DMDALocalInfo &info,
                                 Parameters **P,
                                 const Vector2 ***x,
                                 Vector2 ***R) {

  double
    x_min   = m_grid->x0() - m_grid->Lx(),
    y_min   = m_grid->y0() - m_grid->Ly(),
    dx      = m_grid->dx(),
    dy      = m_grid->dy(),
    scaling = 1.0;

  // Compute the residual at Dirichlet BC nodes and reset the residual to zero elsewhere.
  //
  // here we loop over all the *owned* nodes
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    for (int i = info.xs; i < info.xs + info.xm; i++) {

      // compute the residual at ice-free map plane locations
      if ((int)P[j][i].node_type == NODE_EXTERIOR) {
        for (int k = info.zs; k < info.zs + info.zm; k++) {
          R[j][i][k] = scaling * (x[j][i][k] - u_exterior); // STORAGE_ORDER
        }
        continue;
      }

      for (int k = info.zs; k < info.zs + info.zm; k++) {
        // reset to zero
        R[j][i][k] = 0.0;     // STORAGE_ORDER

        // compute the residual at Dirichlet nodes in verification tests
        if (dirichlet_node(info, {i, j, k})) {
          double
            xx = x_min + i * dx,
            yy = y_min + j * dy,
            b  = P[j][i].bed,
            H  = P[j][i].thickness,
            zz = grid_z(b, H, info.mz, k);
          Vector2 U_bc = u_bc(xx, yy, zz);
          R[j][i][k] = scaling * (x[j][i][k] - U_bc); // STORAGE_ORDER
        }
      }
    }
  }
}

/*!
 * Computes the residual.
 */
void Blatter::compute_residual(DMDALocalInfo *petsc_info,
                               const Vector2 ***X, Vector2 ***R) {
  auto info = grid_transpose(*petsc_info);

  // Stencil width of 1 is not very important, but if info.sw > 1 will lead to more
  // redundant computation (we would be looping over elements that don't contribute to any
  // owned nodes).
  assert(info.sw == 1);

  // horizontal grid spacing is the same on all multigrid levels
  double
    x_min = m_grid->x0() - m_grid->Lx(),
    y_min = m_grid->y0() - m_grid->Ly(),
    dx    = m_grid->dx(),
    dy    = m_grid->dy();

  fem::Q1Element3 element(info, fem::Q13DQuadrature8(),
                          dx, dy, x_min, y_min);

  // Number of nodes per element.
  const int Nk = fem::q13d::n_chi;
  assert(element.n_chi() <= Nk);
  assert(element.n_pts() <= m_Nq);

  // scalar quantities
  double z[Nk];
  double floatation[Nk], sea_level[Nk], bottom_elevation[Nk], ice_thickness[Nk], surface_elevation[Nk];
  double B[Nk], basal_yield_stress[Nk];
  int node_type[Nk];

  // 2D vector quantities
  Vector2 velocity[Nk], R_nodal[Nk];

  // FIXME: this communicates ghosts every time the residual is computed, which is excessive.
  //
  // note: we use m_da below because all multigrid levels use the same 2D grid
  DataAccess<Parameters**> P(m_da, 2, GHOSTED);
  // note: we use info.da below because ice hardness is on the grid corresponding to the
  // current multigrid level
  DataAccess<double***> ice_hardness(info.da, 3, GHOSTED);

  // Compute the residual at Dirichlet nodes and set it to zero elsewhere.
  residual_dirichlet(info, P, X, R);

  // loop over all the elements that have at least one owned node
  for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
    for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {

      // Initialize 2D geometric info at element nodes
      nodal_parameter_values(element, P, i, j,
                             node_type,
                             bottom_elevation,
                             ice_thickness,
                             surface_elevation,
                             sea_level);

      // skip ice-free (exterior) elements
      if (exterior_element(node_type)) {
        continue;
      }

      // loop over elements in a column
      for (int k = info.gzs; k < info.gzs + info.gzm - 1; k++) {

        // Reset element residual to zero in preparation.
        for (int n = 0; n < Nk; ++n) {
          R_nodal[n] = 0.0;
        }

        // Compute z coordinates of the nodes of this element
        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(i, j, k, n);
          z[n] = grid_z(bottom_elevation[n], ice_thickness[n], info.mz, I.k);
        }

        // Compute values of chi, chi_x, chi_y, chi_z and quadrature weights at quadrature
        // points on this physical element
        element.reset(i, j, k, z);

        // Get nodal values of ice velocity.
        {
          element.nodal_values(X, velocity);

          // Take care of Dirichlet BC: don't contribute to Dirichlet nodes and set nodal
          // values of the current iterate to Dirichlet BC values.
          for (int n = 0; n < Nk; ++n) {
            auto I = element.local_to_global(n);
            if (dirichlet_node(info, I)) {
              element.mark_row_invalid(n);
              velocity[n] = u_bc(element.x(n), element.y(n), element.z(n));
            }
          }
        }

        // Get nodal values of ice hardness.
        element.nodal_values((double***)ice_hardness, B);

        // "main" part of the residual
        residual_f(element, velocity, B, R_nodal);

        // the "source term" (driving stress)
        residual_source_term(element, surface_elevation, bottom_elevation, R_nodal);

        // basal boundary
        if (k == 0) {
          for (int n = 0; n < Nk; ++n) {
            auto I = element.local_to_global(n);

            basal_yield_stress[n] = P[I.j][I.i].tauc;
            floatation[n]   = P[I.j][I.i].floatation;
          }

          // use an N*N-point equally-spaced quadrature at grounding lines
          fem::Q1Element3Face *face = grounding_line(floatation) ? &m_face100 : &m_face4;
          // face 4 is the bottom face in fem::q13d::incident_nodes
          face->reset(4, z);

          residual_basal(element, *face, basal_yield_stress, floatation, velocity, R_nodal);
        }

        // lateral boundary
        // loop over all vertical faces (see fem::q13d::incident_nodes for the order)
        for (int f = 0; f < 4; ++f) {
          if (marine_boundary(f, node_type, bottom_elevation, sea_level)) {
            // use an N*N-point equally-spaced quadrature for partially-submerged faces
            fem::Q1Element3Face *face = (partially_submerged_face(f, z, sea_level) ?
                                         &m_face100 : &m_face4);
            face->reset(f, z);

            residual_lateral(element, *face, surface_elevation, z, sea_level, R_nodal);
          }
        } // end of the loop over element faces

        // top boundary (verification tests only)
        if (k == info.mz - 2) {
          // face 5 is the top face in fem::q13d::incident_nodes
          m_face4.reset(5, z);

          residual_surface(element, m_face4, R_nodal);
        }

        element.add_contribution(R_nodal, R);
      } // end of the loop over i
    } // end of the loop over j
  } // end of the loop over k
}

PetscErrorCode Blatter::function_callback(DMDALocalInfo *info,
                                          const Vector2 ***x, Vector2 ***f,
                                          Blatter *solver) {
  try {
    solver->compute_residual(info, x, f);
  } catch (...) {
    MPI_Comm com = solver->grid()->com;
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

} // end of namespace stressbalance
} // end of namespace pism
