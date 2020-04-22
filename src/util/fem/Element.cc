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
#include <cmath>                // std::sqrt()

#include "Element.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace fem {

//! Determinant of a 3x3 matrix
static double det(const double a[3][3]) {
  return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
          a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
          a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
}

//! Cross product of two 3D vectors
static Vector3 cross(const Vector3 &a, const Vector3 &b) {
  return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

// extract a row of a 3x3 matrix
static Vector3 row(const double A[3][3], size_t k) {
  return {A[k][0], A[k][1], A[k][2]};
}

// extract a column of a 3x3 matrix
static Vector3 column(const double A[3][3], size_t k) {
  return {A[0][k], A[1][k], A[2][k]};
}

// dot product of a vector and a [dx, dy, dz] vector in Germ
static double dot(const Vector3 &v, const Germ &a) {
  return a.dx * v.x + a.dy * v.y + a.dz + v.z;
}

//! Invert a 3x3 matrix
static void invert(const double A[3][3], double result[3][3]) {
  const double det_A = det(A);

  assert(det_A != 0.0);

  Vector3
    x0 = column(A, 0),
    x1 = column(A, 1),
    x2 = column(A, 2),
    a  = cross(x1, x2),
    b  = cross(x2, x0),
    c  = cross(x0, x1);

  double A_cofactor[3][3] = {{a.x, a.y, a.z},
                             {b.x, b.y, b.z},
                             {c.x, c.y, c.z}};

  // inverse(A) = 1/det(A) * transpose(A_cofactor)
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // note the transpose on the RHS
      result[i][j] = A_cofactor[j][i] / det_A;
    }
  }
}

//! Compute derivatives with respect to x,y using J^{-1} and derivatives with respect to xi, eta.
static Germ multiply(const double A[3][3], const Germ &v) {
  // FIXME: something is not right here
  return {v.val,
          dot(row(A, 0), v),
          dot(row(A, 1), v),
          dot(row(A, 2), v)};
}

static void set_to_identity(double A[3][3]) {
  A[0][0] = 1.0;
  A[0][1] = 0.0;
  A[0][2] = 0.0;
  A[1][0] = 0.0;
  A[1][1] = 1.0;
  A[1][2] = 0.0;
  A[2][0] = 0.0;
  A[2][1] = 0.0;
  A[2][2] = 1.0;
}

Element::Element(const IceGrid &grid, int Nq, int n_chi, int block_size)
  : m_n_chi(n_chi),
    m_Nq(Nq),
    m_block_size(block_size) {

  // get sub-domain information from the grid:
  auto da = grid.get_dm(1, 0);  // dof = 1, stencil_width = 0
  PetscErrorCode ierr = DMDAGetLocalInfo(*da, &m_grid);
  if (ierr != 0) {
    throw std::runtime_error("Failed to allocate an Element instance");
  }
  // reset da: we don't want to end up depending on it
  m_grid.da = NULL;

  m_germs.resize(Nq * n_chi);
  m_row.resize(block_size);
  m_col.resize(block_size);
}

Element::Element(const DMDALocalInfo &grid_info, int Nq, int n_chi, int block_size)
  : m_n_chi(n_chi),
    m_Nq(Nq),
    m_block_size(block_size) {

  m_grid = grid_info;
  // reset da: we don't want to end up depending on it
  m_grid.da = NULL;

  m_germs.resize(Nq * n_chi);
  m_row.resize(block_size);
  m_col.resize(block_size);
}


Element::~Element() {
  // empty
}

//! Initialize shape function values and quadrature weights of a 2D physical element.
/** Assumes that the Jacobian does not depend on coordinates of the current quadrature point.
 */
void Element::initialize(const double J[3][3],
                         ShapeFunction f,
                         unsigned int n_chi,
                         const std::vector<QuadPoint>& points,
                         const std::vector<double>& W) {

  double J_inv[3][3];
  invert(J, J_inv);

  for (unsigned int q = 0; q < m_Nq; q++) {
    for (unsigned int k = 0; k < n_chi; k++) {
      m_germs[q * m_n_chi + k] = multiply(J_inv, f(k, points[q]));
    }
  }

  m_weights.resize(m_Nq);
  const double J_det = det(J);
  for (unsigned int q = 0; q < m_Nq; q++) {
    m_weights[q] = J_det * W[q];
  }
}

Element2::Element2(const IceGrid &grid, int Nq, int n_chi, int block_size)
  : Element(grid, Nq, n_chi, block_size) {
  // empty
}

Element2::Element2(const DMDALocalInfo &grid_info, int Nq, int n_chi, int block_size)
  : Element(grid_info, Nq, n_chi, block_size) {
  // empty
}

void Element2::nodal_values(const IceModelVec2Int &x_global, int *result) const {
  for (unsigned int k = 0; k < m_n_chi; ++k) {
    const int
      ii = m_i + m_i_offset[k],
      jj = m_j + m_j_offset[k];
    result[k] = x_global.as_int(ii, jj);
  }
}

/*!@brief Initialize the Element to element (`i`, `j`) for the purposes of inserting into
  global residual and Jacobian arrays. */
void Element2::reset(int i, int j) {
  m_i = i;
  m_j = j;

  for (unsigned int k = 0; k < m_n_chi; ++k) {
    m_col[k].i = i + m_i_offset[k];
    m_col[k].j = j + m_j_offset[k];
    m_col[k].k = 0;
    m_col[k].c = 0;
  }
  m_row = m_col;

  // We never sum into rows that are not owned by the local rank.
  for (unsigned int k = 0; k < m_n_chi; k++) {
    int pism_i = 0, pism_j = 0;
    local_to_global(k, pism_i, pism_j);
    if (pism_i < m_grid.xs or m_grid.xs + m_grid.xm - 1 < pism_i or
        pism_j < m_grid.ys or m_grid.ys + m_grid.ym - 1 < pism_j) {
      mark_row_invalid(k);
    }
  }
}

/*!@brief Mark that the row corresponding to local degree of freedom `k` should not be updated
  when inserting into the global residual or Jacobian arrays. */
void Element::mark_row_invalid(int k) {
  m_row[k].i = m_row[k].j = m_row[k].k = m_invalid_dof;
}

/*!@brief Mark that the column corresponding to local degree of freedom `k` should not be updated
  when inserting into the global Jacobian arrays. */
void Element::mark_col_invalid(int k) {
  m_col[k].i = m_col[k].j = m_col[k].k = m_invalid_dof;
}

//! Add the contributions of an element-local Jacobian to the global Jacobian matrix.
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
                                                   m_block_size, m_row.data(),
                                                   m_block_size, m_col.data(),
                                                   K, ADD_VALUES);
  PISM_CHK(ierr, "MatSetValuesBlockedStencil");
}

Q1Element2::Q1Element2(const IceGrid &grid, const Quadrature &quadrature)
  : Element2(grid, quadrature.weights().size(), q1::n_chi, q1::n_chi) {

  double dx = grid.dx();
  double dy = grid.dy();

  m_i_offset = {0, 1, 1, 0};
  m_j_offset = {0, 0, 1, 1};

  // south, east, north, west
  m_normals = {{0.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}};

  m_side_lengths = {dx, dy, dx, dy};

  // compute the Jacobian

  double J[3][3];
  set_to_identity(J);
  J[0][0] = 0.5 * dx;
  J[0][1] = 0.0;
  J[1][0] = 0.0;
  J[1][1] = 0.5 * dy;

  // initialize germs and quadrature weights for the quadrature on this physical element
  initialize(J, q1::chi, q1::n_chi, quadrature.points(), quadrature.weights());
  reset(0, 0);
}

Q1Element2::Q1Element2(const DMDALocalInfo &grid_info,
                       double dx,
                       double dy,
                       const Quadrature &quadrature)
  : Element2(grid_info, quadrature.weights().size(), q1::n_chi, q1::n_chi) {

  m_i_offset = {0, 1, 1, 0};
  m_j_offset = {0, 0, 1, 1};

  // south, east, north, west
  m_normals = {{0.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}};

  m_side_lengths = {dx, dy, dx, dy};

  // compute the Jacobian

  double J[3][3];
  set_to_identity(J);
  J[0][0] = 0.5 * dx;
  J[0][1] = 0.0;
  J[1][0] = 0.0;
  J[1][1] = 0.5 * dy;

  // initialize germs and quadrature weights for the quadrature on this physical element
  initialize(J, q1::chi, q1::n_chi, quadrature.points(), quadrature.weights());
  reset(0, 0);
}

P1Element2::P1Element2(const IceGrid &grid, const Quadrature &quadrature, int type)
  : Element2(grid, quadrature.weights().size(), p1::n_chi, q1::n_chi) {

  double dx = grid.dx();
  double dy = grid.dy();

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
    m_normals = {n12, n20, n01};
    pts = {p[1], p[2], p[0]};
    break;
  case 2:
    m_i_offset = {1, 0, 1};
    m_j_offset = {1, 1, 0};
    m_normals = {n23, -1.0 * n13, n12};
    pts = {p[2], p[3], p[1]};
    break;
  case 3:
  default:
    m_i_offset = {0, 0, 1};
    m_j_offset = {1, 0, 1};
    m_normals = {n30, -1.0 * n20, n23};
    pts = {p[3], p[0], p[2]};
    break;
  }

  reset(0, 0);
  // make sure add_contribution() ignores entries corresponding to the 4-th node
  mark_row_invalid(3);
  mark_col_invalid(3);

  m_side_lengths.resize(n_sides());
  // compute side lengths
  for (unsigned int k = 0; k < n_sides(); ++k) {
    m_side_lengths[k] = (pts[k] - pts[(k + 1) % 3]).magnitude();
  }

  // compute the Jacobian
  double J[3][3];
  set_to_identity(J);
  J[0][0] = pts[1].u - pts[0].u;
  J[0][1] = pts[1].v - pts[0].v;
  J[1][0] = pts[2].u - pts[0].u;
  J[1][1] = pts[2].v - pts[0].v;

  // initialize germs and quadrature weights for the quadrature on this physical element
  initialize(J, p1::chi, p1::n_chi, quadrature.points(), quadrature.weights());
  reset(0, 0);
}

Element3::Element3(const DMDALocalInfo &grid_info, int Nq, int n_chi, int block_size)
  : Element(grid_info, Nq, n_chi, block_size) {
  m_i = 0;
  m_j = 0;
  m_k = 0;
}
Element3::Element3(const IceGrid &grid, int Nq, int n_chi, int block_size)
  : Element(grid, Nq, n_chi, block_size) {
  m_i = 0;
  m_j = 0;
  m_k = 0;
}

Q1Element3::Q1Element3(const DMDALocalInfo &grid_info,
                       double dx,
                       double dy,
                       const Quadrature &quadrature)
  : Element3(grid_info, quadrature.weights().size(), q13d::n_chi, q13d::n_chi),
    m_dx(dx),
    m_dy(dy),
    m_points(quadrature.points()),
    m_w(quadrature.weights()) {

  m_weights.resize(m_Nq);

  m_i_offset = {0, 1, 1, 0, 0, 1, 1, 0};
  m_j_offset = {0, 0, 1, 1, 0, 0, 1, 1};
  m_k_offset = {0, 0, 0, 0, 1, 1, 1, 1};

  // store values of shape functions on the reference element
  m_chi.resize(m_Nq * m_n_chi);
  for (unsigned int q = 0; q < m_Nq; q++) {
    for (unsigned int n = 0; n < m_n_chi; n++) {
      m_chi[q * m_n_chi + n] = q13d::chi(n, m_points[q]);
    }
  }
  m_germs = m_chi;
}

Q1Element3::Q1Element3(const IceGrid &grid, const Quadrature &quadrature)
  : Element3(grid, quadrature.weights().size(), q13d::n_chi, q13d::n_chi),
    m_dx(grid.dx()),
    m_dy(grid.dy()),
    m_points(quadrature.points()),
    m_w(quadrature.weights()) {

  m_weights.resize(m_Nq);

  m_i_offset = {0, 1, 1, 0, 0, 1, 1, 0};
  m_j_offset = {0, 0, 1, 1, 0, 0, 1, 1};
  m_k_offset = {0, 0, 0, 0, 1, 1, 1, 1};

  // store values of shape functions on the reference element
  m_chi.resize(m_Nq * m_n_chi);
  for (unsigned int q = 0; q < m_Nq; q++) {
    for (unsigned int n = 0; n < m_n_chi; n++) {
      m_chi[q * m_n_chi + n] = q13d::chi(n, m_points[q]);
    }
  }
}


/*! Initialize the element `i,j,k`.
 *
 *
 * @param[in] i i-index of the lower left node
 * @param[in] j j-index of the lower left node
 * @param[in] k k-index of the lower left node
 * @param[in] z z-coordinates of the nodes of this element
 */
void Q1Element3::reset(int i, int j, int k, const std::vector<double> &z) {
  // Record i,j,k corresponding to the current element:
  m_i = i;
  m_j = j;
  m_k = k;

  // Set row and column info used to add contributions:
  for (unsigned int n = 0; n < m_n_chi; ++n) {
    m_col[n].i = k + m_k_offset[n]; // x -- z (STORAGE_ORDER)
    m_col[n].j = i + m_i_offset[n]; // y -- x (STORAGE_ORDER)
    m_col[n].k = j + m_j_offset[n]; // z -- y (STORAGE_ORDER)
    m_col[n].c = 0;
  }
  m_row = m_col;

  // Mark rows that we don't own as invalid:
  for (unsigned int n = 0; n < m_n_chi; n++) {
    auto I = local_to_global(n);
    if (I.i < m_grid.xs or m_grid.xs + m_grid.xm - 1 < I.i or
        I.j < m_grid.ys or m_grid.ys + m_grid.ym - 1 < I.j) {
      mark_row_invalid(n);
    }
  }

  // Compute J^{-1} and use it to compute m_germs and m_weights:
  for (unsigned int q = 0; q < m_Nq; q++) {

    Vector3 dz{0.0, 0.0, 0.0};
    for (unsigned int n = 0; n < m_n_chi; ++n) {
      auto &chi = m_chi[q * m_n_chi + n];
      dz.x += chi.dx * z[n];
      dz.y += chi.dy * z[n];
      dz.z += chi.dz * z[n];
    }

    double J[3][3] = {{m_dx / 2.0,        0.0, dz.x},
                      {       0.0, m_dy / 2.0, dz.y},
                      {       0.0,        0.0, dz.z}};

    double J_det = J[0][0] * J[1][1] * J[2][2];

    assert(J_det != 0.0);

    double J_inv[3][3] = {{1.0 / J[0][0],           0.0, -J[0][2] / (J[0][0] * J[2][2])},
                          {0.0,           1.0 / J[1][1], -J[1][2] / (J[1][1] * J[2][2])},
                          {0.0,                     0.0,                 1.0 / J[2][2]}};

    m_weights[q] = J_det * m_w[q];

    for (unsigned int n = 0; n < m_n_chi; n++) {
      auto &chi = m_chi[q * m_n_chi + n];
      // FIXME: I should be able to use multiply() defined above, but there must be a bug
      // there...
      m_germs[q * m_n_chi + n] = {chi.val,
                                  J_inv[0][0] * chi.dx + J_inv[0][1] * chi.dy + J_inv[0][2] * chi.dz,
                                  J_inv[1][0] * chi.dx + J_inv[1][1] * chi.dy + J_inv[1][2] * chi.dz,
                                  J_inv[2][0] * chi.dx + J_inv[2][1] * chi.dy + J_inv[2][2] * chi.dz};
    }
  }
}

Q1Element3Face::Q1Element3Face(double dx, double dy, const Quadrature &quadrature)
  : m_dx(dx),
    m_dy(dy),
    m_points(quadrature.points()),
    m_w(quadrature.weights()),
    m_n_chi(q13d::n_chi),
    m_Nq(m_w.size()) {

  // Note: here I set m_n_chi to q13d::n_chi (8) while each face has only four basis
  // functions that are not zero on it. One could make evaluate() cheaper by omitting
  // these zeros. (This would make the code somewhat more complex.)

  m_chi.resize(m_Nq * m_n_chi);
  m_weights.resize(m_Nq);
  m_normals.resize(m_Nq);
}

void Q1Element3Face::reset(int face, const std::vector<double> &z) {

  // Turn coordinates of a 2D quadrature point on [-1,1]*[-1,1] into coordinates on a face
  // of the cube [-1,1]*[-1,1]*[-1,1].
  for (unsigned int q = 0; q < m_Nq; q++) {
    auto pt = m_points[q];
    QuadPoint P;

    // This expresses parameterizations of faces of the reference element: for example,
    // face 0 is parameterized by (s, t) -> [-1, s, t].
    switch (face) {
    case 0:
      P = {-1.0, pt.xi, pt.eta};
      break;
    case 1:
      P = { 1.0, pt.xi, pt.eta};
      break;
    case 2:
      P = {pt.xi, -1.0, pt.eta};
      break;
    case 3:
      P = {pt.xi,  1.0, pt.eta};
      break;
    case 4:
      P = {pt.xi, pt.eta, -1.0};
      break;
    case 5:
      P = {pt.xi, pt.eta,  1.0};
      break;
    }

    // Compute dz/dxi, dz/deta and dz/dzeta.
    Vector3 dz{0.0, 0.0, 0.0};
    for (unsigned int n = 0; n < m_n_chi; ++n) {
      // Note: chi(n, point) for a particular face does not depend on the element geometry
      // and could be computed in advance (in the constructor). On the other hand, I
      // expect that the number of faces in a Neumann boundary is relatively small and
      // this optimization is not likely to pay off.
      auto chi = q13d::chi(n, P);

      m_chi[q * m_n_chi + n] = chi.val;

      dz.x += chi.dx * z[n];
      dz.y += chi.dy * z[n];
      dz.z += chi.dz * z[n];
    }

    // Use the magnitude of the normal to a face to turn quadrature weights on
    // [-1,1]*[-1,1] into quadrature weights on a face of this physical element.

    m_weights[q] = m_w[q];

    // Sign (-1 or +1) used to set normal orientation This relies of the order of faces in
    // fem::q13d::incident_nodes and above (see the code setting QuadPoint P).
    double sign = 2 * (face % 2) - 1;

    switch (face) {
    case 0:
    case 1:
      m_weights[q] *= 0.5 * m_dy * dz.z;
      m_normals[q] = {sign, 0.0, 0.0};
      break;
    case 2:
    case 3:
      m_weights[q] *= 0.5 * m_dx * dz.z;
      m_normals[q] = {0.0, sign, 0.0};
      break;
    case 4:
    case 5:
      {
        double
          a =  0.5 * m_dy * dz.x,
          b =  0.5 * m_dx * dz.y,
          c = 0.25 * m_dx * m_dy,
          M = std::sqrt(a * a + b * b + c * c);
        m_weights[q] *= M;

        M *= sign;
        m_normals[q] = {a / M, b / M, c / M};
      }
      break;
    }
  } // end of the loop over quadrature points
}


} // end of namespace fem
} // end of namespace pism
