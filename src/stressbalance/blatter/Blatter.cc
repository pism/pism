/* Copyright (C) 2020 PISM Authors
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

#include <cassert>              // assert
#include <cmath>                // std::pow, std::fabs
#include <algorithm>            // std::max
#include <cstring>              // memset

#include "Blatter.hh"
#include "pism/util/fem/FEM.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Vector2.hh"

#include "DataAccess.hh"
#include "grid_hierarchy.hh"
#include "pism/util/node_types.hh"

#include "pism/rheology/FlowLaw.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/basalstrength/basal_resistance.hh"

#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace stressbalance {

const Vector2 u_exterior = {0.0, 0.0};

/*!
 * 2D input parameters
 */
struct Parameters {
  // elevation (z coordinate) of the bottom domain boundary
  double bed;
  // thickness of the domain
  double thickness;
  // NodeType stored as double
  double node_type;
  // basal yield stress
  double tauc;
  // sea level elevation (used to determine if a location is grounded)
  double sea_level;
  // floatation function (positive where floating, zero or negative where grounded)
  double floatation;
};

/*!
 * Compute node type using domain thickness and the thickness threshold `min_thickness`.
 *
 * An element contains ice if ice thickness at all its nodes equal or exceeds the
 * `min_thickness` threshold.
 *
 * A node is *interior* if all four elements it belongs to contain ice.
 *
 * A node is *exterior* if it belongs to zero icy elements.
 *
 * A node that is neither interior nor exterior is a *boundary* node.
 */
static void blatter_node_type(DM da, double min_thickness) {
  // Note that P provides access to a ghosted copy of 2D parameters, so changes to P have
  // no lasting effect.
  DataAccess<Parameters**> P(da, 2, GHOSTED);

  DMDALocalInfo info;
  int ierr = DMDAGetLocalInfo(da, &info); PISM_CHK(ierr, "DMDAGetLocalInfo");
  info = grid_transpose(info);

  // loop over all the owned nodes and reset node type
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    for (int i = info.xs; i < info.xs + info.xm; i++) {
      P[j][i].node_type = 0;
    }
  }

  // Note that dx, dy, and quadrature don't matter here.
  fem::Q1Element2 E(info, 1.0, 1.0, fem::Q1Quadrature1());

  Parameters p[fem::q1::n_chi];

  // Loop over all the elements with at least one owned node and compute the number of icy
  // elements each node belongs to.
  for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
    for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {
      E.reset(i, j);

      E.nodal_values((Parameters**)P, p);

      // An element is "interior" (contains ice) if all of its nodes have thickness above
      // the threshold
      bool interior = true;
      for (int k = 0; k < fem::q1::n_chi; ++k) {
        if (p[k].thickness < min_thickness) {
          interior = false;
          break;
        }
      }

      for (int k = 0; k < fem::q1::n_chi; ++k) {
        int ii, jj;
        E.local_to_global(k, ii, jj);
        P[jj][ii].node_type += interior;
      }
    }
  }

  DataAccess<Parameters**> result(da, 2, NOT_GHOSTED);

  // Loop over all the owned nodes and turn the number of "icy" elements this node belongs
  // to into node type.
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    for (int i = info.xs; i < info.xs + info.xm; i++) {

      switch ((int)P[j][i].node_type) {
      case 4:
        result[j][i].node_type = NODE_INTERIOR;
        break;
      case 0:
        result[j][i].node_type = NODE_EXTERIOR;
        break;
      default:
        result[j][i].node_type = NODE_BOUNDARY;
      }
    }
  }
}

/*!
 * Returns true if a node is in the Dirichlet part of the boundary, false otherwise.
 *
 * Used by verification tests.
 */
static bool dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I) {
  (void) info;
  (void) I;
  return false;
}

/*! Dirichlet BC
*/
static Vector2 u_bc(double x, double y, double z) {
  (void) x;
  (void) y;
  (void) z;
  return {0.0, 0.0};
}

/*!
 * Return true if an element does not contain ice, i.e. is a part of the "exterior" of the
 * ice mass.
 *
 * @param[in] node_type node type at the nodes of an element (an array of 8 integers; only
 *                      4 are used)
 */
static bool exterior_element(const int *node_type) {
  // number of nodes per map-plane cell
  int N = 4;
  for (int n = 0; n < N; ++n) {
    if (node_type[n] == NODE_EXTERIOR) {
      return true;
    }
  }
  return false;
}

/*!
 * Return true if the current map-plane cell contains the grounding line, false otherwise.
 *
 * This is used to determine whether to use more quadrature points to estimate integrals
 * over the bottom face of the basal element.
 *
 * The code takes advantage of the ordering of element nodes: lower 4 first, then upper 4.
 * This means that we can loop over the first 4 nodes and ignore the other 4.
 */
static bool grounding_line(const double *F) {

  // number of nodes per map-plane cell
  int N = 4;

  bool
    grounded = false,
    floating = false;

  for (int n = 0; n < N; ++n) {
    if (F[n] <= 0.0) {
      grounded = true;
    } else {
      floating = true;
    }
  }

  return grounded and floating;
}

/*!
 * Return true if the current vertical face is partially submerged.
 *
 * This is used to determine whether to use more quadrature points to estimate integrals
 * over this face when computing lateral boundary conditions.
 */
static bool partially_submerged_face(int face, const double *z, const double *z_sl) {
  auto nodes = fem::q13d::incident_nodes[face];

  // number of nodes per face
  int N = 4;

  bool
    above = false,
    below = false;

  for (int n = 0; n < N; ++n) {
    int k = nodes[n];
    if (z[k] > z_sl[k]) {
      above = true;
    } else {
      below = true;
    }
  }

  return above and below;
}

/*!
 * Return true if the current face is a part of the Neumann boundary, false otherwise.
 *
 * A face is a part of the Neumann boundary if all four nodes are Neumann nodes. If a node
 * is *both* a Neumann and a Dirichlet node (this may happen), then we treat it as a
 * Neumann node here: element.add_contribution() will do the right thing in this case.
 */
static bool neumann_bc_face(int face, const int *node_type) {
  auto nodes = fem::q13d::incident_nodes[face];

  // number of nodes per face
  int N = 4;
  for (int n = 0; n < N; ++n) {
    if (not (node_type[nodes[n]] == NODE_BOUNDARY)) {
      return false;
    }
  }
  return true;
}

/*! Set the residual at Dirichlet locations
 *
 * Compute the residual at Dirichlet locations and reset the residual to zero elsewhere.
 *
 * Setting it to zero is necessary because we call DMDASNESSetFunctionLocal() with
 * INSERT_VALUES.
 *
 */
static void residual_dirichlet(const GridInfo &grid_info,
                               const DMDALocalInfo &info,
                               Parameters **P,
                               const Vector2 ***x,
                               Vector2 ***R) {
  double
    dx = grid_info.dx(info.mx),
    dy = grid_info.dy(info.my);

  // Compute the residual at Dirichlet BC nodes and reset the residual to zero elsewhere.
  //
  // here we loop over all the *owned* nodes
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    for (int i = info.xs; i < info.xs + info.xm; i++) {
      for (int k = info.zs; k < info.zs + info.zm; k++) {

        // Dirichlet nodes
        if (dirichlet_node(info, {i, j, k}) or
            (int)P[j][i].node_type == NODE_EXTERIOR) {

          // Dirichlet scale
          Vector2 s = {1.0, 1.0};

          Vector2 U_bc;
          if (dirichlet_node(info, {i, j, k})) {
            double
              xx = grid_info.x(dx, i),
              yy = grid_info.y(dy, j),
              b  = P[j][i].bed,
              H  = P[j][i].thickness,
              zz = grid_z(b, H, info.mz, k);
            U_bc = u_bc(xx, yy, zz);
          } else {
            U_bc = u_exterior;
          }

          Vector2 r = x[j][i][k] - U_bc;

          R[j][i][k] = {r.u * s.u, r.v * s.v}; // STORAGE_ORDER
        } else {
          R[j][i][k] = 0.0;     // STORAGE_ORDER
        }
      }
    }
  }
}

void Blatter::compute_residual(DMDALocalInfo *petsc_info,
                               const Vector2 ***x, Vector2 ***R) {
  auto info = grid_transpose(*petsc_info);

  double
    g = m_config->get_number("constants.standard_gravity"),
    rho_w = m_config->get_number("constants.sea_water.density");

  // Stencil width of 1 is not very important, but if info.sw > 1 will lead to more
  // redundant computation (we would be looping over elements that don't contribute to any
  // owned nodes).
  assert(info.sw == 1);

  // Compute grid spacing from domain dimensions and the grid size
  //
  // FIXME: we don't need to re-compute dx and dy now that we don't coarsen in horizontal
  // directions.
  double
    dx = m_grid_info.dx(info.mx),
    dy = m_grid_info.dy(info.my);

  double
    ds_max             = m_config->get_number("stress_balance.blatter.max_cliff_height"),
    grad_s_max         = ds_max / std::max(dx, dy),
    grad_s_squared_max = pow(grad_s_max, 2.0);

  fem::Q1Element3 element(info, dx, dy, fem::Q13DQuadrature8());
  fem::Q1Element3Face
    face4(dx, dy, fem::Q1Quadrature4()),     // 4-point Gaussian quadrature
    face100(dx, dy, fem::Q1QuadratureN(10)); // 100-point quadrature for grounding lines
                                             // and partially-submerged faces

  DataAccess<Parameters**> P(info.da, 2, GHOSTED);
  DataAccess<double***> ice_hardness(info.da, 3, GHOSTED);

  residual_dirichlet(m_grid_info, info, P, x, R);

  // Maximum number of nodes per element.
  const int Nk = fem::q13d::n_chi;
  assert(element.n_chi() <= Nk);

  // Maximum number of quadrature points per element or face
  const int Nq = 100;
  assert(element.n_pts() <= Nq);
  assert(face4.n_pts() <= Nq);
  assert(face100.n_pts() <= Nq);

  // scalar quantities evaluated at quadrature points
  double x_nodal[Nk], xq[Nq];
  double y_nodal[Nk], yq[Nq];
  double B_nodal[Nk], Bq[Nq];
  double s_nodal[Nk], s[Nq], s_x[Nq], s_y[Nq], s_z[Nq];
  double sl_nodal[Nk], z_sl[Nq];
  double tauc_nodal[Nk], tauc[Nq];
  double f_nodal[Nk], floatation[Nq];

  std::vector<double> z_nodal(Nk);
  double zq[Nq];

  // 2D vector quantities evaluated at quadrature points
  Vector2 u_nodal[Nk], u[Nq], u_x[Nq], u_y[Nq], u_z[Nq];

  // quantities evaluated at element nodes
  Vector2 R_nodal[Nk];
  int node_type[Nk];
  double b_nodal[Nk];
  double H_nodal[Nk];

  // loop over all the elements that have at least one owned node
  for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
    for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {

      // Compute or fetch 2D geometric info
      for (int n = 0; n < Nk; ++n) {
        auto I = element.local_to_global(i, j, 0, n);

        auto p = P[I.j][I.i];

        node_type[n] = p.node_type;

        x_nodal[n] = m_grid_info.x(dx, I.i);
        y_nodal[n] = m_grid_info.y(dy, I.j);

        s_nodal[n]  = p.bed + p.thickness;
        sl_nodal[n] = p.sea_level;

        b_nodal[n] = p.bed;
        H_nodal[n] = p.thickness;
      }

      // skip ice-free (exterior) elements
      if (exterior_element(node_type)) {
        continue;
      }

      for (int k = info.gzs; k < info.gzs + info.gzm - 1; k++) {

        // Reset element residual to zero in preparation.
        memset(R_nodal, 0, sizeof(R_nodal));

        // Compute coordinates of the nodes of this element and fetch node types.
        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(i, j, k, n);
          z_nodal[n] = grid_z(b_nodal[n], H_nodal[n], info.mz, I.k);
        }

        // compute values of chi, chi_x, chi_y, chi_z and quadrature weights at quadrature
        // points on this physical element
        element.reset(i, j, k, z_nodal);

        // Get nodal values of ice hardness.
        element.nodal_values((double***)ice_hardness, B_nodal);

        // Get nodal values of u.
        element.nodal_values(x, u_nodal);

        // Take care of Dirichlet BC: don't contribute to Dirichlet nodes and set nodal
        // values of the current iterate to Dirichlet BC values.
        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(n);
          if (dirichlet_node(info, I)) {
            element.mark_row_invalid(n);
            u_nodal[n] = u_bc(x_nodal[n], y_nodal[n], z_nodal[n]);
          }
        }

        // evaluate u and its partial derivatives at quadrature points
        element.evaluate(u_nodal, u, u_x, u_y, u_z);

        // evaluate B (ice hardness) at quadrature points
        element.evaluate(B_nodal, Bq);

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
          m_flow_law->effective_viscosity(Bq[q], gamma, &eta, nullptr);

          // loop over all test functions
          for (int t = 0; t < Nk; ++t) {
            const auto &psi = element.chi(q, t);

            R_nodal[t].u += W * (eta * (psi.dx * (4.0 * ux + 2.0 * vy) +
                                        psi.dy * (uy + vx) +
                                        psi.dz * uz));
            R_nodal[t].v += W * (eta * (psi.dx * (uy + vx) +
                                        psi.dy * (2.0 * ux + 4.0 * vy) +
                                        psi.dz * vz));
          }
        }

        // compute the surface gradient at quadrature points
        element.evaluate(s_nodal, s, s_x, s_y, s_z);

        for (int q = 0; q < element.n_pts(); ++q) {
          auto W = element.weight(q);

          // limit the surface slope to avoid artifacts near steep margins
          Vector2 grad_s(s_x[q], s_y[q]);

          if (grad_s.magnitude_squared() > grad_s_squared_max) {
            grad_s = (grad_s_max / grad_s.magnitude()) * grad_s;
          }

          // loop over all test functions
          for (int t = 0; t < Nk; ++t) {
            const auto &psi = element.chi(q, t);

            R_nodal[t].u += W * psi.val * m_rhog * grad_s.u;
            R_nodal[t].v += W * psi.val * m_rhog * grad_s.v;
          }
        }

        // include basal drag
        if (k == 0) {

          for (int n = 0; n < Nk; ++n) {
            auto I = element.local_to_global(n);

            tauc_nodal[n] = P[I.j][I.i].tauc;
            f_nodal[n]    = P[I.j][I.i].floatation;
          }

          // use an N*N-point equally-spaced quadrature at grounding lines
          fem::Q1Element3Face *face = grounding_line(f_nodal) ? &face100 : &face4;

          // face 4 is the bottom face in fem::q13d::incident_nodes
          face->reset(4, z_nodal);

          face->evaluate(u_nodal, u);
          face->evaluate(tauc_nodal, tauc);
          face->evaluate(f_nodal, floatation);

          for (int q = 0; q < face->n_pts(); ++q) {
            auto W = face->weight(q);

            bool grounded = floatation[q] <= 0.0;
            double beta = grounded ? m_basal_sliding_law->drag(tauc[q], u[q].u, u[q].v) : 0.0;

            // loop over all test functions
            for (int t = 0; t < Nk; ++t) {
              auto psi = face->chi(q, t);

              R_nodal[t].u += W * psi * beta * u[q].u;
              R_nodal[t].v += W * psi * beta * u[q].v;
            }
          }
        }

        // loop over all vertical faces (see fem::q13d::incident_nodes for the order)
        for (int f = 0; f < 4; ++f) {

          if (neumann_bc_face(f, node_type)) {
            // use an N*N-point equally-spaced quadrature at for partially-submerged faces
            fem::Q1Element3Face *face = partially_submerged_face(f, z_nodal.data(), sl_nodal) ? &face100 : &face4;

            face->reset(f, z_nodal);

            // compute physical coordinates of quadrature points on this face
            face->evaluate(x_nodal, xq);
            face->evaluate(y_nodal, yq);
            face->evaluate(z_nodal.data(), zq);
            face->evaluate(sl_nodal, z_sl);

            // loop over all quadrature points
            for (int q = 0; q < face->n_pts(); ++q) {
              auto W = face->weight(q);
              auto N3 = face->normal(q);
              Vector2 N = {N3.x, N3.y};

              double p = rho_w * g * std::max(z_sl[q] - zq[q], 0.0);

              // loop over all test functions
              for (int t = 0; t < Nk; ++t) {
                auto psi = face->chi(q, t);

                R_nodal[t] += W * psi * p * N;
              }
            }
          } // end of "if (neumann_bc_face())"
        } // end of the loop over element faces

        element.add_contribution(R_nodal, R);
      } // end of the loop over i
    } // end of the loop over j
  } // end of the loop over k
}

/*!
 * Set the Jacobian to identity at Dirichlet nodes.
 */
static void jacobian_dirichlet(const DMDALocalInfo &info, Parameters **P, Mat J) {
  PetscErrorCode ierr;

  // take care of Dirichlet nodes (both explicit and grid points outside the domain)
  //
  // here we loop over all the *owned* nodes
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    for (int i = info.xs; i < info.xs + info.xm; i++) {
      for (int k = info.zs; k < info.zs + info.zm; k++) {
        if ((int)P[j][i].node_type == NODE_EXTERIOR or dirichlet_node(info, {i, j, k})) {

          // Dirichlet scaling
          Vector2 scaling = {1.0, 1.0};
          double identity[4] = {scaling.u, 0, 0, scaling.v};

          MatStencil row;
          row.i = k;            // STORAGE_ORDER
          row.j = i;            // STORAGE_ORDER
          row.k = j;            // STORAGE_ORDER
          ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, identity, ADD_VALUES);
          PISM_CHK(ierr, "MatSetValuesBlockedStencil"); // this may throw
        }
      }
    }
  }

}

void Blatter::compute_jacobian(DMDALocalInfo *petsc_info,
                               const Vector2 ***x, Mat A, Mat J) {
  auto info = grid_transpose(*petsc_info);

  (void) x;

  // Zero out the Jacobian in preparation for updating it.
  PetscErrorCode ierr = MatZeroEntries(J);
  PISM_CHK(ierr, "MatZeroEntries");

  // Stencil width of 1 is not very important, but if info.sw > 1 will lead to more
  // redundant computation (we would be looping over elements that don't contribute to any
  // owned nodes).
  assert(info.sw == 1);

  // Compute grid spacing from domain dimensions and the grid size
  //
  // FIXME: we don't need to re-compute dx and dy now that we don't coarsen in horizontal
  // directions.
  double
    dx = m_grid_info.dx(info.mx),
    dy = m_grid_info.dy(info.my);

  fem::Q1Element3 element(info, dx, dy, fem::Q13DQuadrature8());

  fem::Q1Element3Face
    face4(dx, dy, fem::Q1Quadrature4()),     // 4-point Gaussian quadrature
    face100(dx, dy, fem::Q1QuadratureN(10)); // 100-point quadrature for grounding lines

  DataAccess<Parameters**> P(info.da, 2, GHOSTED);
  DataAccess<double***> hardness(info.da, 3, GHOSTED);

  // Maximum number of nodes per element
  const int Nk = fem::q13d::n_chi;
  assert(element.n_chi() <= Nk);

  // Maximum number of quadrature points per element or face
  const int Nq = 100;
  assert(element.n_pts() <= Nq);
  assert(face4.n_pts() <= Nq);
  assert(face100.n_pts() <= Nq);

  // 2D vector quantities evaluated at quadrature points
  Vector2 u_nodal[Nk], u[Nq], u_x[Nq], u_y[Nq], u_z[Nq];

  // scalar quantities evaluated at quadrature points
  double B_nodal[Nk], Bq[Nq];
  double tauc_nodal[Nk], tauc[Nq];
  double f_nodal[Nk], floatation[Nq];

  // scalar quantities evaluated at element nodes
  int node_type[Nk];
  double x_nodal[Nk];
  double y_nodal[Nk];
  std::vector<double> z_nodal(Nk);

  // loop over all the elements that have at least one owned node
  for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
    for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {
      for (int k = info.gzs; k < info.gzs + info.gzm - 1; k++) {

        // Element-local Jacobian matrix (there are Nk vector valued degrees of freedom
        // per element, for a total of Nk*Nk = 64 entries in the local Jacobian.
        double K[2*Nk][2*Nk];
        memset(K, 0, sizeof(K));

        // Compute coordinates of the nodes of this element and fetch node types.
        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(i, j, k, n);

          auto p = P[I.j][I.i];

          node_type[n] = p.node_type;

          x_nodal[n] = m_grid_info.x(dx, I.i);
          y_nodal[n] = m_grid_info.y(dy, I.j);
          z_nodal[n] = grid_z(p.bed, p.thickness, info.mz, I.k);
        }

        // skip ice-free (exterior) elements
        if (exterior_element(node_type)) {
          continue;
        }

        // compute values of chi, chi_x, chi_y, chi_z and quadrature weights at quadrature
        // points on this physical element
        element.reset(i, j, k, z_nodal);

        // Get nodal values of u.
        element.nodal_values(x, u_nodal);

        // Don't contribute to Dirichlet nodes
        for (int n = 0; n < Nk; ++n) {
          auto I = element.local_to_global(n);
          if (dirichlet_node(info, I)) {
            element.mark_row_invalid(n);
            element.mark_col_invalid(n);
            u_nodal[n] = u_bc(x_nodal[n], y_nodal[n], z_nodal[n]);
          }
        }

        // evaluate partial derivatives at quadrature points
        element.evaluate(u_nodal, u, u_x, u_y, u_z);

        // evaluate hardness at quadrature points
        element.nodal_values((double***)hardness, B_nodal);
        element.evaluate(B_nodal, Bq);

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

          double eta, deta;
          m_flow_law->effective_viscosity(Bq[q], gamma, &eta, &deta);

          // loop over test and trial functions, computing the upper-triangular part of
          // the element Jacobian
          for (int t = 0; t < Nk; ++t) {
            auto psi = element.chi(q, t);
            for (int s = t; s < Nk; ++s) {
              auto phi = element.chi(q, s);

              double
                gamma_u = 2.0 * ux * phi.dx + vy * phi.dx + 0.5 * phi.dy * (uy + vx) + 0.5 * uz * phi.dz,
                gamma_v = 2.0 * vy * phi.dy + ux * phi.dy + 0.5 * phi.dx * (uy + vx) + 0.5 * vz * phi.dz;

              double
                eta_u = deta * gamma_u,
                eta_v = deta * gamma_v;

              // Picard part
              K[t * 2 + 0][s * 2 + 0] += W * eta * (4.0 * psi.dx * phi.dx + psi.dy * phi.dy + psi.dz * phi.dz);
              K[t * 2 + 0][s * 2 + 1] += W * eta * (2.0 * psi.dx * phi.dy + psi.dy * phi.dx);
              K[t * 2 + 1][s * 2 + 0] += W * eta * (2.0 * psi.dy * phi.dx + psi.dx * phi.dy);
              K[t * 2 + 1][s * 2 + 1] += W * eta * (4.0 * psi.dy * phi.dy + psi.dx * phi.dx + psi.dz * phi.dz);
              // extra Newton terms
              K[t * 2 + 0][s * 2 + 0] += W * eta_u * (psi.dx * (4.0 * ux + 2.0 * vy) +
                                                      psi.dy * (uy + vx) +
                                                      psi.dz * uz);
              K[t * 2 + 0][s * 2 + 1] += W * eta_v * (psi.dx * (4.0 * ux + 2.0 * vy) +
                                                      psi.dy * (uy + vx) +
                                                      psi.dz * uz);
              K[t * 2 + 1][s * 2 + 0] += W * eta_u * (psi.dx * (uy + vx) +
                                                      psi.dy * (4.0 * vy + 2.0 * ux) +
                                                      psi.dz * vz);
              K[t * 2 + 1][s * 2 + 1] += W * eta_v * (psi.dx * (uy + vx) +
                                                      psi.dy * (4.0 * vy + 2.0 * ux) +
                                                      psi.dz * vz);
            }
          }
        } // end of the loop over q

        // include basal drag
        if (k == 0) {

          for (int n = 0; n < Nk; ++n) {
            auto I = element.local_to_global(n);

            tauc_nodal[n] = P[I.j][I.i].tauc;
            f_nodal[n]    = P[I.j][I.i].floatation;
          }

          fem::Q1Element3Face *face = grounding_line(f_nodal) ? &face100 : &face4;

          // face 4 is the bottom face in fem::q13d::incident_nodes
          face->reset(4, z_nodal);

          face->evaluate(u_nodal, u);
          face->evaluate(tauc_nodal, tauc);
          face->evaluate(f_nodal, floatation);

          for (int q = 0; q < face->n_pts(); ++q) {
            auto W = face->weight(q);

            bool grounded = floatation[q] <= 0.0;
            double beta = 0.0, dbeta = 0.0;
            if (grounded) {
              m_basal_sliding_law->drag_with_derivative(tauc[q], u[q].u, u[q].v, &beta, &dbeta);
            }

            // loop over all test functions
            for (int t = 0; t < Nk; ++t) {
              auto psi = face->chi(q, t);
              for (int s = 0; s < Nk; ++s) {
                auto phi = face->chi(q, s);

                double p = psi * phi;

                K[t * 2 + 0][s * 2 + 0] += W * p * (beta + dbeta * u[q].u * u[q].u);
                K[t * 2 + 0][s * 2 + 1] += W * p * dbeta * u[q].u * u[q].v;
                K[t * 2 + 1][s * 2 + 0] += W * p * dbeta * u[q].v * u[q].u;
                K[t * 2 + 1][s * 2 + 1] += W * p * (beta + dbeta * u[q].v * u[q].v);
              }
            }
          }
        }

        // fill the lower-triangular part of the element Jacobian using the fact that J is
        // symmetric
        for (int t = 0; t < Nk; ++t) {
          for (int s = 0; s < t; ++s) {
            K[t * 2 + 0][s * 2 + 0] = K[s * 2 + 0][t * 2 + 0];
            K[t * 2 + 1][s * 2 + 0] = K[s * 2 + 0][t * 2 + 1];
            K[t * 2 + 0][s * 2 + 1] = K[s * 2 + 1][t * 2 + 0];
            K[t * 2 + 1][s * 2 + 1] = K[s * 2 + 1][t * 2 + 1];
          }
        }

        element.add_contribution(&K[0][0], J);
      } // end of the loop over i
    } // end of the loop over j
  } // end of the loop over k

  jacobian_dirichlet(info, P, J);

  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); PISM_CHK(ierr, "MatAssemblyBegin");
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); PISM_CHK(ierr, "MatAssemblyEnd");
  if (A != J) {
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); PISM_CHK(ierr, "MatAssemblyBegin");
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); PISM_CHK(ierr, "MatAssemblyEnd");
  }

  ierr = MatSetOption(J, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");

  ierr = MatSetOption(J, MAT_SYMMETRIC, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");
}

PetscErrorCode Blatter::function_callback(DMDALocalInfo *info,
                                           const Vector2 ***x, Vector2 ***f,
                                           CallbackData *data) {
  try {
    data->solver->compute_residual(info, x, f);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)data->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

PetscErrorCode Blatter::jacobian_callback(DMDALocalInfo *info,
                                           const Vector2 ***x,
                                           Mat A, Mat J, CallbackData *data) {
  try {
    data->solver->compute_jacobian(info, x, A, J);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)data->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

/*!
 * Allocate the Blatter-Pattyn stress balance solver
 *
 * @param[in] grid PISM's grid.
 * @param[in] Mz number of vertical levels
 * @param[in] n_levels maximum number of grid levels to use
 * @param[in] coarsening_factor grid coarsening factor
 */
Blatter::Blatter(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor)
  : ShallowStressBalance(grid) {

  auto pism_da = grid->get_dm(1, 0);

  int ierr = setup(*pism_da, Mz, n_levels, coarsening_factor);
  if (ierr != 0) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Failed to allocate a Blatter solver instance");
  }

  {
    int mz = Mz + grid_padding(Mz, coarsening_factor, n_levels);
    std::vector<double> sigma(mz);
    double dz = 1.0 / (mz - 1);
    for (int i = 0; i < mz; ++i) {
      sigma[i] = i * dz;
    }
    sigma.back() = 1.0;

    std::map<std::string,std::string> z_attrs =
      {{"axis", "Z"},
       {"long_name", "scaled Z-coordinate in the ice (z_base=0, z_surface=1)"},
       {"units", "1"},
       {"positive", "up"}};

    m_u_sigma.reset(new IceModelVec3(grid, "uvel_sigma", "z_sigma", sigma, z_attrs));
    m_u_sigma->set_attrs("diagnostic", "u velocity component on the sigma grid", "m s-1", "m s-1", "", 0);

    m_v_sigma.reset(new IceModelVec3(grid, "vvel_sigma", "z_sigma", sigma, z_attrs));
    m_v_sigma->set_attrs("diagnostic", "v velocity component on the sigma grid", "m s-1", "m s-1", "", 0);
  }

  {
    rheology::FlowLawFactory ice_factory("stress_balance.blatter.", m_config, m_EC);
    ice_factory.remove(ICE_GOLDSBY_KOHLSTEDT);
    m_flow_law = ice_factory.create();
  }

  m_rhog = (m_config->get_number("constants.ice.density") *
            m_config->get_number("constants.standard_gravity"));
}

/*!
 * Restrict 2D and 3D model parameters from a fine grid to a coarse grid.
 *
 * Re-compute node types from geometry.
 *
 * This hook is called every time SNES needs to update coarse-grid data.
 *
 * FIXME: parameters restricted by this hook do not change from one SNES iteration to the
 * next, so we can return early after the first one.
 */
static PetscErrorCode blatter_restriction_hook(DM fine,
                                               Mat mrestrict, Vec rscale, Mat inject,
                                               DM coarse, void *ctx) {
  // Get rid of warnings about unused arguments
  (void) mrestrict;
  (void) rscale;
  (void) inject;
  GridInfo *grid_info = (GridInfo*)ctx;

  PetscErrorCode ierr;

  // FIXME: we don't need to restrict in 2D now that we don't coarsen in horizontal
  // directions.
  ierr = restrict_data(fine, coarse, "2D_DM"); CHKERRQ(ierr);

  ierr = restrict_data(fine, coarse, "3D_DM"); CHKERRQ(ierr);

  blatter_node_type(coarse, grid_info->min_thickness);

  return 0;
}

PetscErrorCode blatter_coarsening_hook(DM dm_fine, DM dm_coarse, void *ctx) {
  PetscErrorCode ierr;

  ierr = setup_level(dm_coarse, *(GridInfo*)ctx); CHKERRQ(ierr);

  ierr = DMCoarsenHookAdd(dm_coarse, blatter_coarsening_hook, blatter_restriction_hook, ctx); CHKERRQ(ierr);

  // 2D
  //
  // FIXME: we don't need the 2D restriction now that we don't coarsen in horizontal
  // directions.
  ierr = create_restriction(dm_fine, dm_coarse, "2D_DM"); CHKERRQ(ierr);

  // 3D
  ierr = create_restriction(dm_fine, dm_coarse, "3D_DM"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode Blatter::setup(DM pism_da, int Mz, int n_levels, int coarsening_factor) {
  PetscErrorCode ierr;
  // DM
  //
  // Note: in the PISM's DA pism_da PETSc's and PISM's meaning of x and y are the same.
  {
    PetscInt dim, Mx, My, Nx, Ny;
    PetscInt
      Nz            = 1,
      dof           = 2,        // u and v velocity components
      stencil_width = 1;

    ierr = DMDAGetInfo(pism_da,
                       &dim,
                       &Mx,
                       &My,
                       NULL,             // Mz
                       &Nx,              // number of processors in y-direction
                       &Ny,              // number of processors in x-direction
                       NULL,             // ditto, z-direction
                       NULL,             // number of degrees of freedom per node
                       NULL,             // stencil width
                       NULL, NULL, NULL, // types of ghost nodes at the boundary
                       NULL);            // stencil width
    CHKERRQ(ierr);

    assert(dim == 2);

    const PetscInt *lx, *ly;
    ierr = DMDAGetOwnershipRanges(pism_da, &lx, &ly, NULL); CHKERRQ(ierr);

    double
      x_max = m_grid->Lx(),
      x_min = -x_max,
      y_max = m_grid->Ly(),
      y_min = -y_max;

    // pad the vertical grid to allow for n_levels multigrid levels
    Mz += grid_padding(Mz, coarsening_factor, n_levels);

    ierr = DMDACreate3d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, // STORAGE_ORDER
                        DMDA_STENCIL_BOX,
                        Mz, Mx, My,                         // STORAGE_ORDER
                        Nz, Nx, Ny,                         // STORAGE_ORDER
                        dof,                                // dof
                        stencil_width,                      // stencil width
                        NULL, lx, ly,                       // STORAGE_ORDER
                        m_da.rawptr()); CHKERRQ(ierr);

    // semi-coarsening: coarsen in the vertical direction only
    ierr = DMDASetRefinementFactor(m_da, coarsening_factor, 1, 1); CHKERRQ(ierr); // STORAGE_ORDER

    ierr = DMSetFromOptions(m_da); CHKERRQ(ierr);

    ierr = DMSetUp(m_da); CHKERRQ(ierr);

    double min_thickness = m_config->get_number("stress_balance.ice_free_thickness_standard");

    m_grid_info = {x_min, x_max,
                   y_min, y_max,
                   min_thickness,
                   sizeof(Parameters)/sizeof(double)};

    // set up 2D and 3D parameter storage
    ierr = setup_level(m_da, m_grid_info); CHKERRQ(ierr);

    // tell PETSc how to coarsen this grid and how to restrict data to a coarser grid
    ierr = DMCoarsenHookAdd(m_da, blatter_coarsening_hook, blatter_restriction_hook, &m_grid_info);
    CHKERRQ(ierr);
  }

  // Vec
  {
    ierr = DMCreateGlobalVector(m_da, m_x.rawptr()); CHKERRQ(ierr);
  }

  // SNES
  {
    ierr = SNESCreate(m_grid->com, m_snes.rawptr()); CHKERRQ(ierr);

    // ierr = SNESSetOptionsPrefix(m_snes, "blatter_"); CHKERRQ(ierr);

    ierr = SNESSetDM(m_snes, m_da); CHKERRQ(ierr);

    m_callback_data.da = m_da;
    m_callback_data.solver = this;

    ierr = DMDASNESSetFunctionLocal(m_da, INSERT_VALUES,
                                    (DMDASNESFunction)function_callback,
                                    &m_callback_data); CHKERRQ(ierr);

    ierr = DMDASNESSetJacobianLocal(m_da,
                                    (DMDASNESJacobian)jacobian_callback,
                                    &m_callback_data); CHKERRQ(ierr);

    ierr = SNESSetFromOptions(m_snes); CHKERRQ(ierr);
  }

  return 0;
}

/*!
 * Set 2D parameters on the finest grid.
 */
void Blatter::init_2d_parameters(const Inputs &inputs) {

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    alpha         = ice_density / water_density;

  const IceModelVec2S
    &tauc      = *inputs.basal_yield_stress,
    &H         = inputs.geometry->ice_thickness,
    &b         = inputs.geometry->bed_elevation,
    &sea_level = inputs.geometry->sea_level_elevation;

  {
    DataAccess<Parameters**> P(m_da, 2, NOT_GHOSTED);

    IceModelVec::AccessList list{&tauc, &H, &b, &sea_level};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double
        b_grounded = b(i, j),
        b_floating = sea_level(i, j) - alpha * H(i, j),
        s_grounded = b(i, j) + H(i, j),
        s_floating = sea_level(i, j) + (1.0 - alpha) * H(i, j);

      P[j][i].tauc = tauc(i, j);
      P[j][i].thickness = H(i, j);
      P[j][i].sea_level = sea_level(i, j);
      P[j][i].bed = std::max(b_grounded, b_floating);
      P[j][i].node_type = NODE_EXTERIOR;
      P[j][i].floatation = s_floating - s_grounded;
    }
  }

  blatter_node_type(m_da, m_grid_info.min_thickness);
}

/*!
 * Set 3D parameters on the finest grid.
 */
void Blatter::init_ice_hardness(const Inputs &inputs) {

  auto enthalpy = inputs.enthalpy;

  DMDALocalInfo info;
  int ierr = DMDAGetLocalInfo(m_da, &info); PISM_CHK(ierr, "DMDAGetLocalInfo");
  info = grid_transpose(info);

  const auto &zlevels = enthalpy->levels();
  auto Mz = zlevels.size();

  DataAccess<Parameters**> P2(m_da, 2, NOT_GHOSTED);
  DataAccess<double***> P3(m_da, 3, NOT_GHOSTED);

  IceModelVec::AccessList list{enthalpy};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double H  = P2[j][i].thickness;

    const double *E = enthalpy->get_column(i, j);

    for (int k = 0; k < info.mz; ++k) {
      double
        z        = grid_z(0.0, H, info.mz, k),
        depth    = H - z,
        pressure = m_EC->pressure(depth),
        E_local  = 0.0;

      auto k0 = m_grid->kBelowHeight(z);

      if (k0 + 1 < Mz) {
        double lambda = (z - zlevels[k0]) / (zlevels[k0 + 1] - zlevels[k0]);

        E_local = (1.0 - lambda) * E[k0] + lambda * E[k0 + 1];
      } else {
        E_local = E[Mz - 1];
      }

      P3[j][i][k] = m_flow_law->hardness(E_local, pressure);
    }

  } // end of the loop over grid points
}

Blatter::~Blatter() {
  // empty
}

void Blatter::init_impl() {
  m_log->message(2, "* Initializing the Blatter stress balance...\n");

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {
    File input_file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
    bool u_sigma_found = input_file.find_variable("uvel_sigma");
    bool v_sigma_found = input_file.find_variable("vvel_sigma");
    unsigned int start = input_file.nrecords() - 1;

    if (u_sigma_found and v_sigma_found) {
      m_log->message(3, "Reading uvel_sigma and vvel_sigma...\n");

      m_u_sigma->read(input_file, start);
      m_v_sigma->read(input_file, start);

      set_initial_guess(*m_u_sigma, *m_v_sigma);
    } else {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "uvel_sigma and vvel_sigma not found");
    }
  } else {
    int ierr = VecSet(m_x, 0.0); PISM_CHK(ierr, "VecSet");
  }
}

void Blatter::define_model_state_impl(const File &output) const {
  m_u_sigma->define(output);
  m_v_sigma->define(output);
}

void Blatter::write_model_state_impl(const File &output) const {
  m_u_sigma->write(output);
  m_v_sigma->write(output);
}

void Blatter::update(const Inputs &inputs, bool full_update) {
  (void) inputs;
  (void) full_update;

  init_2d_parameters(inputs);
  init_ice_hardness(inputs);

  int ierr = SNESSolve(m_snes, NULL, m_x); PISM_CHK(ierr, "SNESSolve");

  // report the number of iterations
  {
    PetscInt            its, lits;
    SNESConvergedReason reason;
    SNESGetIterationNumber(m_snes, &its);
    SNESGetConvergedReason(m_snes, &reason);
    SNESGetLinearSolveIterations(m_snes, &lits);
    m_log->message(2, "%s: SNES: %d, KSP: %d\n",
                   SNESConvergedReasons[reason], (int)its, (int)lits);
  }
  // FIXME: check if SNESSolve() succeeded and try to recover

  // put basal velocity in m_velocity to use it in the next call
  get_basal_velocity(m_velocity);
  compute_basal_frictional_heating(m_velocity, *inputs.basal_yield_stress,
                                   inputs.geometry->cell_type,
                                   m_basal_frictional_heating);

  compute_averaged_velocity(m_velocity);

  // copy the solution from m_x to m_u_sigma, m_v_sigma for re-starting
  copy_solution();
}

void Blatter::copy_solution() {
  Vector2 ***x = nullptr;
  int ierr = DMDAVecGetArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecGetArray");

  int Mz = m_u_sigma->levels().size();

  IceModelVec::AccessList list{m_u_sigma.get(), m_v_sigma.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto u = m_u_sigma->get_column(i, j);
    auto v = m_v_sigma->get_column(i, j);

    for (int k = 0; k < Mz; ++k) {
      u[k] = x[j][i][k].u;      // STORAGE_ORDER
      v[k] = x[j][i][k].v;      // STORAGE_ORDER
    }
  }

  ierr = DMDAVecRestoreArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecRestoreArray");
}

void Blatter::get_basal_velocity(IceModelVec2V &result) {
  Vector2 ***x = nullptr;
  int ierr = DMDAVecGetArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecGetArray");

  IceModelVec::AccessList list{&result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j).u = x[j][i][0].u;      // STORAGE_ORDER
    result(i, j).v = x[j][i][0].v;      // STORAGE_ORDER
  }

  ierr = DMDAVecRestoreArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecRestoreArray");
}


void Blatter::set_initial_guess(const IceModelVec3 &u_sigma,
                                const IceModelVec3 &v_sigma) {
  Vector2 ***x = nullptr;
  int ierr = DMDAVecGetArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecGetArray");

  int Mz = m_u_sigma->levels().size();

  IceModelVec::AccessList list{&u_sigma, &v_sigma};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto u = u_sigma.get_column(i, j);
    auto v = v_sigma.get_column(i, j);

    for (int k = 0; k < Mz; ++k) {
      x[j][i][k].u = u[k];      // STORAGE_ORDER
      x[j][i][k].v = v[k];      // STORAGE_ORDER
    }
  }

  ierr = DMDAVecRestoreArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecRestoreArray");
}

void Blatter::compute_averaged_velocity(IceModelVec2V &result) {
  PetscErrorCode ierr;

  Vector2 ***x = nullptr;
  ierr = DMDAVecGetArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecGetArray");

  int Mz = m_u_sigma->levels().size();

  IceModelVec::AccessList list{&result};
  DataAccess<Parameters**> P2(m_da, 2, NOT_GHOSTED);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double H = P2[j][i].thickness;

    Vector2 V(0.0, 0.0);

    if (H > 0.0) {
      // use trapezoid rule to compute the column average
      double dz = H / (Mz - 1);
      for (int k = 0; k < Mz - 1; ++k) {
        V += x[j][i][k] + x[j][i][k + 1]; // STORAGE_ORDER
      }
      V *= (0.5 * dz) / H;
    }

    result(i, j) = V;
  }

  ierr = DMDAVecRestoreArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecRestoreArray");

  result.update_ghosts();
}


IceModelVec3::Ptr Blatter::velocity_u_sigma() const {
  return m_u_sigma;
}

IceModelVec3::Ptr Blatter::velocity_v_sigma() const {
  return m_v_sigma;
}

} // end of namespace stressbalance
} // end of namespace pism
