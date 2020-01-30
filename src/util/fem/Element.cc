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
#include "Element.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace fem {

//! Determinant of a square matrix of size 2.
static double determinant(const double J[2][2]) {
  return J[0][0] * J[1][1] - J[1][0] * J[0][1];
}

//! Compute the inverse of a two by two matrix.
static void invert(const double A[2][2], double A_inv[2][2]) {
  const double det_A = determinant(A);

  assert(det_A != 0.0);

  A_inv[0][0] =  A[1][1] / det_A;
  A_inv[0][1] = -A[0][1] / det_A;
  A_inv[1][0] = -A[1][0] / det_A;
  A_inv[1][1] =  A[0][0] / det_A;
}

//! Compute derivatives with respect to x,y using J^{-1} and derivatives with respect to xi, eta.
static Germ multiply(const double A[2][2], const Germ &v) {
  Germ result;
  result.val = v.val;
  result.dx  = v.dx * A[0][0] + v.dy * A[0][1];
  result.dy  = v.dx * A[1][0] + v.dy * A[1][1];
  return result;
}

// Multiply a matrix by a vector.
static Vector2 multiply(const double A[2][2], const Vector2 &v) {
  return {v.u * A[0][0] + v.v * A[0][1],
          v.u * A[1][0] + v.v * A[1][1]};
}

Element::Element(const IceGrid &grid, const Quadrature &quadrature)
  : m_n_chi(0),
    m_Nq(quadrature.points().size()),
    m_grid(grid) {

  m_germs = (Germs*) malloc(m_Nq * m_n_chi_max * sizeof(Germ));
  if (not m_germs) {
    throw std::runtime_error("Failed to allocate an Element instance");
  }

  reset(0, 0);
}

Element::~Element() {
  free(m_germs);
  m_germs = nullptr;
}

//! Initialize shape function values and quadrature weights of a 2D physical element.
/** Assumes that the Jacobian does not depend on coordinates of the current quadrature point.
 */
void Element::initialize(ShapeFunction2 f,
                         unsigned int n_chi,
                         const std::vector<QuadPoint>& points,
                         const std::vector<double>& W) {

  double J_inv[2][2];
  invert(m_J, J_inv);

  for (unsigned int q = 0; q < m_Nq; q++) {
    for (unsigned int k = 0; k < n_chi; k++) {
      m_germs[q][k] = multiply(J_inv, f(k, points[q]));
    }
  }

  m_weights.resize(m_Nq);
  const double J_det = determinant(m_J);
  for (unsigned int q = 0; q < m_Nq; q++) {
    m_weights[q] = J_det * W[q];
  }
}

void Element::nodal_values(const IceModelVec2Int &x_global, int *result) const {
  for (unsigned int k = 0; k < m_n_chi; ++k) {
    const int
      ii = m_i + m_i_offset[k],
      jj = m_j + m_j_offset[k];
    result[k] = x_global.as_int(ii, jj);
  }
}

/*!@brief Initialize the Element to element (`i`, `j`) for the purposes of inserting into
  global residual and Jacobian arrays. */
void Element::reset(int i, int j) {
  m_i = i;
  m_j = j;

  for (unsigned int k = 0; k < m_n_chi; ++k) {
    m_col[k].i = i + m_i_offset[k];
    m_col[k].j = j + m_j_offset[k];
    m_col[k].k = 0;

    m_row[k].i = m_col[k].i;
    m_row[k].j = m_col[k].j;
    m_row[k].k = m_col[k].k;
  }

  // We never sum into rows that are not owned by the local rank.
  for (unsigned int k = 0; k < m_n_chi; k++) {
    int pism_i = m_row[k].i, pism_j = m_row[k].j;
    if (pism_i < m_grid.xs() or m_grid.xs() + m_grid.xm() - 1 < pism_i or
        pism_j < m_grid.ys() or m_grid.ys() + m_grid.ym() - 1 < pism_j) {
      mark_row_invalid(k);
    }
  }
}

/*!@brief Mark that the row corresponding to local degree of freedom `k` should not be updated
  when inserting into the global residual or Jacobian arrays. */
void Element::mark_row_invalid(int k) {
  m_row[k].i = m_row[k].j = m_invalid_dof;
}

/*!@brief Mark that the column corresponding to local degree of freedom `k` should not be updated
  when inserting into the global Jacobian arrays. */
void Element::mark_col_invalid(int k) {
  m_col[k].i = m_col[k].j = m_invalid_dof;
}

//! Add the contributions of an element-local Jacobian to the global Jacobian vector.
/*! The element-local Jacobian should be given as a row-major array of
 *  Nk*Nk values in the scalar case or (2Nk)*(2Nk) values in the
 *  vector valued case.
 *
 *  Note that MatSetValuesBlockedStencil ignores negative indexes, so
 *  values in K corresponding to locations marked using
 *  mark_row_invalid() and mark_col_invalid() are ignored. (Just as they
 *  should be.)
 */
void Element::add_contribution(const double *K, Mat J) const {
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,
                                                   m_n_chi, m_row,
                                                   m_n_chi, m_col,
                                                   K, ADD_VALUES);
  PISM_CHK(ierr, "MatSetValuesBlockedStencil");
}

Q1Element::Q1Element(const IceGrid &grid, const Quadrature &quadrature)
  : Element(grid, quadrature) {

  double dx = grid.dx();
  double dy = grid.dy();

  m_i_offset = {0, 1, 1, 0};
  m_j_offset = {0, 0, 1, 1};
  m_n_chi = q1::n_chi;

  // south, east, north, west
  m_normals = {{0.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}};

  m_side_lengths = {dx, dy, dx, dy};

  // compute the Jacobian
  m_J[0][0] = 0.5 * dx;
  m_J[0][1] = 0.0;
  m_J[1][0] = 0.0;
  m_J[1][1] = 0.5 * dy;

  // initialize germs and quadrature weights for the quadrature on this physical element
  initialize(q1::chi, q1::n_chi, quadrature.points(), quadrature.weights());
}

P1Element::P1Element(const IceGrid &grid, const Quadrature &quadrature, int type)
  : Element(grid, quadrature) {

  double dx = grid.dx();
  double dy = grid.dy();

  m_n_chi = p1::n_chi;

  // outward pointing normals for all sides of a Q1 element with sides aligned with X and
  // Y axes
  const Vector2
    n01( 0.0, -1.0),  // south
    n12( 1.0,  0.0),  // east
    n23( 0.0,  1.0),  // north
    n30(-1.0,  0.0);  // west

  // "diagonal" sides
  Vector2
    n13( 1.0, dx / dy), // 1-3 diagonal, outward for element 0
    n20(-1.0, dx / dy); // 2-0 diagonal, outward for element 1

  // normalize
  n13 /= n13.magnitude();
  n20 /= n20.magnitude();

  // coordinates of nodes of a Q1 element this P1 element is embedded in (up to
  // translation)
  Vector2 p[] = {{0, 0}, {dx, 0}, {dx, dy}, {0, dy}};
  std::vector<Vector2> pts;

  switch (type) {
  case 0:
    m_i_offset = {0, 1, 0};
    m_j_offset = {0, 0, 1};
    m_normals = {n01, n13, n30};
    pts = {p[0], p[1], p[3]};
    break;
  case 1:
    m_i_offset = {1, 1, 0};
    m_j_offset = {0, 1, 0};
    m_normals = {n01, n12, n20};
    pts = {p[1], p[2], p[0]};
    break;
  case 2:
    m_i_offset = {1, 0, 1};
    m_j_offset = {1, 1, 0};
    m_normals = {n12, n23, -1.0 * n13};
    pts = {p[2], p[3], p[0]};
    break;
  case 3:
  default:
    m_i_offset = {0, 0, 1};
    m_j_offset = {1, 0, 1};
    m_normals = {n23, n30, -1.0 * n20};
    pts = {p[3], p[0], p[2]};
    break;
  }

  m_side_lengths.resize(n_sides());
  // compute side lengths
  for (unsigned int k = 0; k < n_sides(); ++k) {
    m_side_lengths[k] = (pts[k] - pts[(k + 1) % 3]).magnitude();
  }

  // compute the Jacobian
  m_J[0][0] = pts[1].u - pts[0].u;
  m_J[0][1] = pts[1].v - pts[0].v;
  m_J[1][0] = pts[2].u - pts[0].u;
  m_J[1][1] = pts[2].v - pts[0].v;

  // initialize germs and quadrature weights for the quadrature on this physical element
  initialize(p1::chi, p1::n_chi, quadrature.points(), quadrature.weights());
}

} // end of namespace fem
} // end of namespace pism
