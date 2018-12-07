// Copyright (C) 2009--2011, 2013, 2014, 2015, 2016, 2017, 2018 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

// Utility functions used by the SSAFEM code.

#include <cassert>              // assert
#include <cstring>              // memset
#include <cstdlib>              // malloc, free

#include "FETools.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Logger.hh"

#include "pism/util/error_handling.hh"

namespace pism {

//! FEM (Finite Element Method) utilities
namespace fem {

namespace q0 {
/*!
 * Piecewise-constant shape functions.
 */
Germ chi(unsigned int k, const QuadPoint &pt) {
  assert(k < q0::n_chi);

  Germ result;

  if ((k == 0 and pt.xi <= 0.0 and pt.eta <= 0.0) or
      (k == 1 and pt.xi > 0.0 and pt.eta <= 0.0) or
      (k == 2 and pt.xi > 0.0 and pt.eta > 0.0) or
      (k == 3 and pt.xi <= 0.0 and pt.eta > 0.0)) {
    result.val = 1.0;
  } else {
    result.val = 0.0;
  }

  result.dx  = 0.0;
  result.dy  = 0.0;

  return result;
}
} // end of namespace q0

namespace q1 {

//! Q1 basis functions on the reference element with nodes (-1,-1), (1,-1), (1,1), (-1,1).
Germ chi(unsigned int k, const QuadPoint &pt) {
  assert(k < q1::n_chi);

  const double xi[] = {-1.0,  1.0,  1.0, -1.0};
  const double eta[] = {-1.0, -1.0,  1.0,  1.0};

  Germ result;

  result.val = 0.25 * (1.0 + xi[k] * pt.xi) * (1.0 + eta[k] * pt.eta);
  result.dx  = 0.25 * xi[k] * (1.0 + eta[k] * pt.eta);
  result.dy  = 0.25 * eta[k] * (1.0 + xi[k] * pt.xi);

  return result;
}


const Vector2* outward_normals() {
  static const Vector2 n[] = {Vector2( 0.0, -1.0),  // south
                              Vector2( 1.0,  0.0),  // east
                              Vector2( 0.0,  1.0),  // north
                              Vector2(-1.0,  0.0)}; // west
  return n;
}

} // end of namespace q1

namespace p1 {

//! P1 basis functions on the reference element with nodes (0,0), (1,0), (0,1).
Germ chi(unsigned int k, const QuadPoint &pt) {
  assert(k < q1::n_chi);
  Germ result;

  switch (k) {
  case 0:
    result.val = 1.0 - pt.xi - pt.eta;
    result.dx  = -1.0;
    result.dy  = -1.0;
    break;
  case 1:
    result.val = pt.xi;
    result.dx  = 1.0;
    result.dy  = 0.0;
    break;
  case 2:
    result.val = pt.eta;
    result.dx  = 0.0;
    result.dy  = 1.0;
    break;
  default:                      // the fourth (dummy) basis function
    result.val = 0.0;
    result.dx  = 0.0;
    result.dy  = 0.0;
  }
  return result;
}
} // end of namespace p1


ElementIterator::ElementIterator(const IceGrid &g) {
  // Start by assuming ghost elements exist in all directions.
  // Elements are indexed by their lower left vertex.  If there is a ghost
  // element on the right, its i-index will be the same as the maximum
  // i-index of a non-ghost vertex in the local grid.
  xs = g.xs() - 1;                    // Start at ghost to the left.
  int xf = g.xs() + g.xm() - 1; // End at ghost to the right.
  ys = g.ys() - 1;                    // Start at ghost at the bottom.
  int yf = g.ys() + g.ym() - 1; // End at ghost at the top.

  lxs = g.xs();
  int lxf = lxs + g.xm() - 1;
  lys = g.ys();
  int lyf = lys + g.ym() - 1;

  // Now correct if needed. The only way there will not be ghosts is if the
  // grid is not periodic and we are up against the grid boundary.

  if (!(g.periodicity() & X_PERIODIC)) {
    // Leftmost element has x-index 0.
    if (xs < 0) {
      xs = 0;
    }
    // Rightmost vertex has index g.Mx-1, so the rightmost element has index g.Mx-2
    if (xf > (int)g.Mx() - 2) {
      xf  = g.Mx() - 2;
      lxf = g.Mx() - 2;
    }
  }

  if (!(g.periodicity() & Y_PERIODIC)) {
    // Bottom element has y-index 0.
    if (ys < 0) {
      ys = 0;
    }
    // Topmost vertex has index g.My - 1, so the topmost element has index g.My - 2
    if (yf > (int)g.My() - 2) {
      yf  = g.My() - 2;
      lyf = g.My() - 2;
    }
  }

  // Tally up the number of elements in each direction
  xm  = xf - xs + 1;
  ym  = yf - ys + 1;
  lxm = lxf - lxs + 1;
  lym = lyf - lys + 1;
}

ElementMap::ElementMap(const IceGrid &grid)
  : m_grid(grid) {
  reset(0, 0);
}

ElementMap::~ElementMap() {
  // empty
}

void ElementMap::nodal_values(const IceModelVec2Int &x_global, int *result) const {
  for (unsigned int k = 0; k < q1::n_chi; ++k) {
    const int
      ii = m_i + m_i_offset[k],
      jj = m_j + m_j_offset[k];
    result[k] = x_global.as_int(ii, jj);
  }
}

/*!@brief Initialize the ElementMap to element (`i`, `j`) for the purposes of inserting into
  global residual and Jacobian arrays. */
void ElementMap::reset(int i, int j) {
  m_i = i;
  m_j = j;

  for (unsigned int k = 0; k < fem::q1::n_chi; ++k) {
    m_col[k].i = i + m_i_offset[k];
    m_col[k].j = j + m_j_offset[k];
    m_col[k].k = 0;

    m_row[k].i = m_col[k].i;
    m_row[k].j = m_col[k].j;
    m_row[k].k = m_col[k].k;
  }

  // We do not ever sum into rows that are not owned by the local rank.
  for (unsigned int k = 0; k < fem::q1::n_chi; k++) {
    int pism_i = m_row[k].i, pism_j = m_row[k].j;
    if (pism_i < m_grid.xs() || m_grid.xs() + m_grid.xm() - 1 < pism_i ||
        pism_j < m_grid.ys() || m_grid.ys() + m_grid.ym() - 1 < pism_j) {
      mark_row_invalid(k);
    }
  }
}

/*!@brief Mark that the row corresponding to local degree of freedom `k` should not be updated
  when inserting into the global residual or Jacobian arrays. */
void ElementMap::mark_row_invalid(int k) {
  m_row[k].i = m_row[k].j = m_invalid_dof;
  // We are solving a 2D system, so MatStencil::k is not used. Here we
  // use it to mark invalid rows.
  m_row[k].k = 1;
}

/*!@brief Mark that the column corresponding to local degree of freedom `k` should not be updated
  when inserting into the global Jacobian arrays. */
void ElementMap::mark_col_invalid(int k) {
  m_col[k].i = m_col[k].j = m_invalid_dof;
  // We are solving a 2D system, so MatStencil::k is not used. Here we
  // use it to mark invalid columns.
  m_col[k].k = 1;
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
void ElementMap::add_contribution(const double *K, Mat J) const {
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,
                                                   fem::q1::n_chi, m_row,
                                                   fem::q1::n_chi, m_col,
                                                   K, ADD_VALUES);
  PISM_CHK(ierr, "MatSetValuesBlockedStencil");
}

const int ElementMap::m_i_offset[4] = {0, 1, 1, 0};
const int ElementMap::m_j_offset[4] = {0, 0, 1, 1};

Quadrature::Quadrature(unsigned int N)
  : m_Nq(N) {

  m_W = (double*) malloc(m_Nq * sizeof(double));
  if (m_W == NULL) {
    throw std::runtime_error("Failed to allocate a Quadrature instance");
  }

  m_germs = (Germs*) malloc(m_Nq * q1::n_chi * sizeof(Germ));
  if (m_germs == NULL) {
    free(m_W);
    throw std::runtime_error("Failed to allocate a Quadrature instance");
  }
}

Germ Quadrature::test_function_values(unsigned int q, unsigned int k) const {
  return m_germs[q][k];
}

double Quadrature::weights(unsigned int q) const {
  return m_W[q];
}

UniformQxQuadrature::UniformQxQuadrature(unsigned int size, double dx, double dy, double scaling)
  : Quadrature(size) {
  // We use uniform Cartesian coordinates, so the Jacobian is constant and diagonal on every
  // element.
  //
  // Note that the reference element is [-1,1]^2, hence the extra factor of 1/2.
  m_J[0][0] = 0.5 * dx / scaling;
  m_J[0][1] = 0.0;
  m_J[1][0] = 0.0;
  m_J[1][1] = 0.5 * dy / scaling;
}

Quadrature::~Quadrature() {
  free(m_W);
  m_W = NULL;

  free(m_germs);
  m_germs = NULL;
}

//! Build quadrature points and weights for a tensor product quadrature based on a 1D quadrature
//! rule. Uses the same 1D quadrature in both directions.
/**
   @param[in] n 1D quadrature size (the resulting quadrature has size n*n)
   @param[in] points1 1D quadrature points
   @param[in] weights1 1D quadrature weights
   @param[out] points resulting 2D quadrature points
   @param[out] weights resulting 2D quadrature weights
 */
static void tensor_product_quadrature(unsigned int n,
                                      const double *points1,
                                      const double *weights1,
                                      QuadPoint *points,
                                      double *weights) {
  unsigned int q = 0;
  for (unsigned int j = 0; j < n; ++j) {
    for (unsigned int i = 0; i < n; ++i) {
      points[q].xi = points1[i];
      points[q].eta = points1[j];

      weights[q] = weights1[i] * weights1[j];

      ++q;
    }
  }
}

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
static Germ apply_jacobian_inverse(const double J_inv[2][2], const Germ &f) {
  Germ result;
  result.val = f.val;
  result.dx  = f.dx * J_inv[0][0] + f.dy * J_inv[0][1];
  result.dy  = f.dx * J_inv[1][0] + f.dy * J_inv[1][1];
  return result;
}

//! Two-by-two Gaussian quadrature on a rectangle.
Q1Quadrature4::Q1Quadrature4(double dx, double dy, double L)
  : UniformQxQuadrature(m_size, dx, dy, L) {

  // coordinates and weights of the 2-point 1D Gaussian quadrature
  const double
    A           = 1.0 / sqrt(3.0),
    points2[2]  = {-A, A},
    weights2[2] = {1.0, 1.0};

  QuadPoint points[m_size];
  double W[m_size];

  tensor_product_quadrature(2, points2, weights2, points, W);

  initialize(q1::chi, q1::n_chi, points, W);
}

Q1Quadrature9::Q1Quadrature9(double dx, double dy, double L)
  : UniformQxQuadrature(m_size, dx, dy, L) {
  // The quadrature points on the reference square.

  const double
    A         = 0.0,
    B         = sqrt(0.6),
    points3[3] = {-B, A, B};

  const double
    w1         = 5.0 / 9.0,
    w2         = 8.0 / 9.0,
    weights3[3] = {w1, w2, w1};

  QuadPoint points[m_size];
  double W[m_size];

  tensor_product_quadrature(3, points3, weights3, points, W);

  initialize(q1::chi, q1::n_chi, points, W);
}

Q1Quadrature16::Q1Quadrature16(double dx, double dy, double L)
  : UniformQxQuadrature(m_size, dx, dy, L) {
  // The quadrature points on the reference square.
  const double
    A          = sqrt(3.0 / 7.0 - (2.0 / 7.0) * sqrt(6.0 / 5.0)), // smaller magnitude
    B          = sqrt(3.0 / 7.0 + (2.0 / 7.0) * sqrt(6.0 / 5.0)), // larger magnitude
    points4[4] = {-B, -A, A, B};

  // The weights w_i for Gaussian quadrature on the reference element with these quadrature points
  const double
    w1          = (18.0 + sqrt(30.0)) / 36.0, // larger
    w2          = (18.0 - sqrt(30.0)) / 36.0, // smaller
    weights4[4] = {w2, w1, w1, w2};

  QuadPoint points[m_size];
  double W[m_size];

  tensor_product_quadrature(4, points4, weights4, points, W);

  initialize(q1::chi, q1::n_chi, points, W);
}


//! @brief 1e4-point (100x100) uniform (*not* Gaussian) quadrature for integrating discontinuous
//! functions.
Q0Quadrature1e4::Q0Quadrature1e4(double dx, double dy, double L)
  : UniformQxQuadrature(m_size, dx, dy, L) {
  assert(m_size_1d * m_size_1d == m_size);

  double xi[m_size_1d], w[m_size_1d];
  const double dxi = 2.0 / m_size_1d;
  for (unsigned int k = 0; k < m_size_1d; ++k) {
    xi[k] = -1.0 + dxi*(k + 0.5);
    w[k]  = 2.0 / m_size_1d;
  }

  QuadPoint points[m_size];
  double W[m_size];

  tensor_product_quadrature(m_size_1d, xi, w, points, W);

  initialize(q0::chi, q0::n_chi, points, W);
}

//! @brief 1e4-point (100x100) uniform (*not* Gaussian) quadrature for integrating discontinuous
//! functions.
Q1Quadrature1e4::Q1Quadrature1e4(double dx, double dy, double L)
  : UniformQxQuadrature(m_size, dx, dy, L) {
  assert(m_size_1d * m_size_1d == m_size);

  double xi[m_size_1d], w[m_size_1d];
  const double dxi = 2.0 / m_size_1d;
  for (unsigned int k = 0; k < m_size_1d; ++k) {
    xi[k] = -1.0 + dxi*(k + 0.5);
    w[k]  = 2.0 / m_size_1d;
  }

  QuadPoint points[m_size];
  double W[m_size];

  tensor_product_quadrature(m_size_1d, xi, w, points, W);

  initialize(q1::chi, q1::n_chi, points, W);
}

//! Initialize shape function values and weights of a 2D quadrature.
/** Assumes that the Jacobian does not depend on coordinates of the current quadrature point.
 */
void Quadrature::initialize(ShapeFunction2 f,
                            unsigned int n_chi,
                            const QuadPoint *points,
                            const double *W) {

  double J_inv[2][2];
  invert(m_J, J_inv);

  for (unsigned int q = 0; q < m_Nq; q++) {
    for (unsigned int k = 0; k < n_chi; k++) {
      m_germs[q][k] = apply_jacobian_inverse(J_inv, f(k, points[q]));
    }
  }

  const double J_det = determinant(m_J);
  for (unsigned int q = 0; q < m_Nq; q++) {
    m_W[q] = J_det * W[q];
  }
}

/** Create a quadrature on a P1 element aligned with coordinate axes and embedded in a Q1 element.

    There are four possible P1 elements in a Q1 element. The argument `N` specifies which one,
    numbering them by the node at the right angle in the "reference" element (0,0) -- (1,0) --
    (0,1).
 */
P1Quadrature::P1Quadrature(unsigned int size, unsigned int N,
                           double dx, double dy, double L)
  : Quadrature(size) {

  // Compute the Jacobian. The nodes of the selected triangle are
  // numbered, the unused node is marked with an "X". In all triangles
  // nodes are numbered in the counter-clockwise direction.
  switch (N) {
  case 0:
    /*
    2------X
    |      |
    |      |
    0------1
    */
    m_J[0][0] = dx / L;
    m_J[0][1] = 0.0;
    m_J[1][0] = 0.0;
    m_J[1][1] = dy / L;
    break;
  case 1:
    /*
    X------1
    |      |
    |      |
    2------0
    */
    m_J[0][0] = 0.0;
    m_J[0][1] = dy / L;
    m_J[1][0] = -dx / L;
    m_J[1][1] = 0.0;
    break;
  case 2:
    /*
    1------0
    |      |
    |      |
    X------2
    */
    m_J[0][0] = -dx / L;
    m_J[0][1] = 0.0;
    m_J[1][0] = 0.0;
    m_J[1][1] = -dy / L;
    break;
  case 3:
    /*
    0------2
    |      |
    |      |
    1------X
    */
    m_J[0][0] = 0.0;
    m_J[0][1] = -dy / L;
    m_J[1][0] = dx / L;
    m_J[1][1] = 0.0;
    break;
  }
}

//! Permute shape functions stored in `f` *in place*.
static void permute(const unsigned int p[q1::n_chi], Germ f[q1::n_chi]) {
  Germ tmp[q1::n_chi];

  // store permuted values to avoid overwriting f too soon
  for (unsigned int k = 0; k < q1::n_chi; ++k) {
    tmp[k] = f[p[k]];
  }

  // copy back into f
  for (unsigned int k = 0; k < q1::n_chi; ++k) {
    f[k] = tmp[k];
  }
}

P1Quadrature3::P1Quadrature3(unsigned int N,
                             double dx, double dy, double L)
  : P1Quadrature(m_size, N, dx, dy, L) {

  const double
    one_over_six   = 1.0 / 6.0,
    two_over_three = 2.0 / 3.0;

  const QuadPoint points[3] = {{two_over_three, one_over_six},
                               {one_over_six,   two_over_three},
                               {one_over_six,   one_over_six}};

  const double W[3] = {one_over_six, one_over_six, one_over_six};

  // Note that we use q1::n_chi here.
  initialize(p1::chi, q1::n_chi, points, W);

  // Permute shape function values according to N, the index of this triangle in the Q1 element.
  const unsigned int
    X       = 3,                // index of the dummy shape function
    p[4][4] = {{0, 1, X, 2},
               {2, 0, 1, X},
               {X, 2, 0, 1},
               {1, X, 2, 0}};
  for (unsigned int q = 0; q < m_size; ++q) {
    permute(p[N], m_germs[q]);
  }
}

DirichletData::DirichletData()
  : m_indices(NULL), m_weight(1.0) {
  for (unsigned int k = 0; k < q1::n_chi; ++k) {
    m_indices_e[k] = 0;
  }
}

DirichletData::~DirichletData() {
  finish(NULL);
}

void DirichletData::init(const IceModelVec2Int *indices,
                         const IceModelVec *values,
                         double weight) {
  m_weight = weight;

  if (indices != NULL) {
    indices->begin_access();
    m_indices = indices;
  }

  if (values != NULL) {
    values->begin_access();
  }
}

void DirichletData::finish(const IceModelVec *values) {
  if (m_indices != NULL) {
    MPI_Comm com = m_indices->grid()->ctx()->com();
    try {
      m_indices->end_access();
      m_indices = NULL;
    } catch (...) {
      handle_fatal_errors(com);
    }
  }

  if (values != NULL) {
    MPI_Comm com = values->grid()->ctx()->com();
    try {
      values->end_access();
    } catch (...) {
      handle_fatal_errors(com);
    }
  }
}

//! @brief Constrain `element`, i.e. ensure that quadratures do not contribute to Dirichlet nodes by marking corresponding rows and columns as "invalid".
void DirichletData::constrain(ElementMap &element) {
  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < q1::n_chi; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      // Mark any kind of Dirichlet node as not to be touched
      element.mark_row_invalid(k);
      element.mark_col_invalid(k);
    }
  }
}

// Scalar version

DirichletData_Scalar::DirichletData_Scalar(const IceModelVec2Int *indices,
                                           const IceModelVec2S *values,
                                           double weight)
  : m_values(values) {
  init(indices, m_values, weight);
}

void DirichletData_Scalar::enforce(const ElementMap &element, double* x_nodal) {
  assert(m_values != NULL);

  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < q1::n_chi; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i = 0, j = 0;
      element.local_to_global(k, i, j);
      x_nodal[k] = (*m_values)(i, j);
    }
  }
}

void DirichletData_Scalar::enforce_homogeneous(const ElementMap &element, double* x_nodal) {
  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < q1::n_chi; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_nodal[k] = 0.;
    }
  }
}

void DirichletData_Scalar::fix_residual(double const *const *const x_global, double **r_global) {
  assert(m_values != NULL);

  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i] = m_weight * (x_global[j][i] - (*m_values)(i,j));
    }
  }
}

void DirichletData_Scalar::fix_residual_homogeneous(double **r_global) {
  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i] = 0.0;
    }
  }
}

void DirichletData_Scalar::fix_jacobian(Mat J) {
  const IceGrid &grid = *m_indices->grid();

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.

  const double identity = m_weight;
  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*m_indices)(i, j) > 0.5) {
        MatStencil row;
        row.j = j;
        row.i = i;
        PetscErrorCode ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, &identity,
                                                         ADD_VALUES);
        PISM_CHK(ierr, "MatSetValuesBlockedStencil"); // this may throw
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

DirichletData_Scalar::~DirichletData_Scalar() {
  finish(m_values);
  m_values = NULL;
}

// Vector version

DirichletData_Vector::DirichletData_Vector(const IceModelVec2Int *indices,
                                           const IceModelVec2V *values,
                                           double weight)
  : m_values(values) {
  init(indices, m_values, weight);
}

void DirichletData_Vector::enforce(const ElementMap &element, Vector2* x_nodal) {
  assert(m_values != NULL);

  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < q1::n_chi; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i = 0, j = 0;
      element.local_to_global(k, i, j);
      x_nodal[k] = (*m_values)(i, j);
    }
  }
}

void DirichletData_Vector::enforce_homogeneous(const ElementMap &element, Vector2* x_nodal) {
  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < q1::n_chi; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_nodal[k].u = 0.0;
      x_nodal[k].v = 0.0;
    }
  }
}

void DirichletData_Vector::fix_residual(Vector2 const *const *const x_global, Vector2 **r_global) {
  assert(m_values != NULL);

  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i] = m_weight * (x_global[j][i] - (*m_values)(i, j));
    }
  }
}

void DirichletData_Vector::fix_residual_homogeneous(Vector2 **r_global) {
  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i].u = 0.0;
      r_global[j][i].v = 0.0;
    }
  }
}

void DirichletData_Vector::fix_jacobian(Mat J) {
  const IceGrid &grid = *m_indices->grid();

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.

  const double identity[4] = {m_weight, 0,
                              0, m_weight};
  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*m_indices)(i, j) > 0.5) {
        MatStencil row;
        row.j = j;
        row.i = i;
        PetscErrorCode ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, identity,
                                                         ADD_VALUES);
        PISM_CHK(ierr, "MatSetValuesBlockedStencil"); // this may throw
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

DirichletData_Vector::~DirichletData_Vector() {
  finish(m_values);
  m_values = NULL;
}

BoundaryQuadrature2::BoundaryQuadrature2(double dx, double dy, double L) {

  const double J[2][2] = {{0.5 * dx / L, 0.0},
                          {0.0, 0.5 * dy / L}};

  // The inverse of the Jacobian.
  double J_inv[2][2];
  invert(J, J_inv);

  // Note that all quadrature weights are 1.0 (and so they are implicitly included below).
  //
  // bottom
  m_W[0] = J[0][0];
  // right
  m_W[1] = J[1][1];
  // top
  m_W[2] = J[0][0];
  // left
  m_W[3] = J[1][1];

  const double C = 1.0 / sqrt(3);
  const QuadPoint points[q1::n_sides][m_Nq] = {
    {{  -C, -1.0}, {   C, -1.0}}, // South
    {{ 1.0,   -C}, { 1.0,    C}}, // East
    {{  -C,  1.0}, {   C,  1.0}}, // North
    {{-1.0,   -C}, {-1.0,    C}}  // West
  };

  memset(m_germs, 0, q1::n_sides*m_Nq*q1::n_chi*sizeof(Germ));

  for (unsigned int side = 0; side < q1::n_sides; ++side) {
    for (unsigned int q = 0; q < m_Nq; ++q) {
      for (unsigned int k = 0; k < q1::n_chi; ++k) {
        Germ phi = q1::chi(k, points[side][q]);

        m_germs[side][q][k] = apply_jacobian_inverse(J_inv, phi);
      }
    }
  }
}

} // end of namespace fem
} // end of namespace pism
