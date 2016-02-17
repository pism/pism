// Copyright (C) 2009--2011, 2013, 2014, 2015, 2016 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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
#include <cassert>
#include <cstring>

#include "FETools.hh"
#include "base/util/IceGrid.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/Logger.hh"

namespace pism {

//! FEM (Finite Element Method) utilities
namespace fem {

const double ShapeQ1::m_xi[]  = {-1.0,  1.0,  1.0, -1.0};
const double ShapeQ1::m_eta[] = {-1.0, -1.0,  1.0,  1.0};

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

void ElementMap::nodal_values(const IceModelVec2Int &x_global, int *result) const
{
  for (unsigned int k = 0; k < ShapeQ1::Nk; ++k) {
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

  for (unsigned int k = 0; k < fem::ShapeQ1::Nk; ++k) {
    m_col[k].i = i + m_i_offset[k];
    m_col[k].j = j + m_j_offset[k];
    m_col[k].k = 0;

    m_row[k].i = m_col[k].i;
    m_row[k].j = m_col[k].j;
    m_row[k].k = m_col[k].k;
  }

  // We do not ever sum into rows that are not owned by the local rank.
  for (unsigned int k = 0; k < fem::ShapeQ1::Nk; k++) {
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
void ElementMap::add_jacobian_contribution(const double *K, Mat J) const {
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,
                                                   fem::ShapeQ1::Nk, m_row,
                                                   fem::ShapeQ1::Nk, m_col,
                                                   K, ADD_VALUES);
  PISM_CHK(ierr, "MatSetValuesBlockedStencil");
}

const int ElementMap::m_i_offset[4] = {0, 1, 1, 0};
const int ElementMap::m_j_offset[4] = {0, 0, 1, 1};

Quadrature_Scalar::Quadrature_Scalar(double dx, double dy, double L)
  : Quadrature2x2(dx, dy, L) {
}

//! Obtain the weights @f$ w_q @f$ for quadrature.
const double* Quadrature2x2::weighted_jacobian() const {
  return m_JxW;
}

//! Obtain the weights @f$ w_q @f$ for quadrature.
Quadrature2x2::Quadrature2x2(double dx, double dy, double L) {
  // Since we use uniform cartesian coordinates, the Jacobian is
  // constant and diagonal on every element.
  //
  // Note that the reference element is @f$ [-1,1]^2 @f$ hence the
  // extra factor of 1/2.
  double jacobian_x = 0.5*dx / L;
  double jacobian_y = 0.5*dy / L;
  m_jacobianDet = jacobian_x*jacobian_y;

  for (unsigned int q = 0; q < Nq; q++) {
    for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
      m_germs[q][k] = ShapeQ1::eval(k, quadPoints[q][0], quadPoints[q][1]);
      m_germs[q][k].dx /= jacobian_x;
      m_germs[q][k].dy /= jacobian_y;
    }
  }

  for (unsigned int q = 0; q < Nq; q++) {
    m_JxW[q] = m_jacobianDet * quadWeights[q];
  }
}

Quadrature_Vector::Quadrature_Vector(double dx, double dy, double L)
  : Quadrature2x2(dx, dy, L) {
}

//! Return the values at all quadrature points of all shape functions.
//* The return value is an Nq by Nk array of Germs. */
const Quadrature2x2::Germs* Quadrature2x2::test_function_values() const {
  return m_germs;
}

/*! @brief Compute the values at the quadrature ponits of a scalar-valued
  finite-element function with element-local degrees of freedom `x_nodal`.*/
/*! There should be room for Quadrature2x2::Nq values in the output vector `vals`. */
void Quadrature_Scalar::quadrature_point_values(const double *x_nodal, double *vals) {
  for (unsigned int q = 0; q < Nq; q++) {
    const Germ *test = m_germs[q];
    vals[q] = 0;
    for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
      vals[q] += test[k].val * x_nodal[k];
    }
  }
}

/*! @brief Compute the values and first derivatives at the quadrature
  points of a scalar-valued finite-element function with element-local
  degrees of freedom `x_nodal`.*/
/*! There should be room for Quadrature2x2::Nq values in the output vectors `vals`, `dx`,
  and `dy`. */
void Quadrature_Scalar::quadrature_point_values(const double *x_nodal,
                                                double *vals, double *dx, double *dy) {
  for (unsigned int q = 0; q < Nq; q++) {
    const Germ *test = m_germs[q];
    vals[q] = 0;
    dx[q] = 0;
    dy[q] = 0;
    for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
      vals[q] += test[k].val * x_nodal[k];
      dx[q]   += test[k].dx * x_nodal[k];
      dy[q]   += test[k].dy * x_nodal[k];
    }
  }
}

/*! @brief Compute the values at the quadrature points of a vector-valued
  finite-element function with element-local degrees of freedom `x_nodal`.*/
/*! There should be room for Quadrature2x2::Nq values in the output vector `vals`. */
void Quadrature_Vector::quadrature_point_values(const Vector2 *x_nodal, Vector2 *result) {
  for (unsigned int q = 0; q < Nq; q++) {
    result[q].u = 0;
    result[q].v = 0;
    const Germ *test = m_germs[q];
    for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
      result[q].u += test[k].val * x_nodal[k].u;
      result[q].v += test[k].val * x_nodal[k].v;
    }
  }
}

/*! @brief Compute the values and symmetric gradient at the quadrature
 *         points of a vector-valued finite-element function with
 *         element-local degrees of freedom `x_nodal`.
 *
 * There should be room for Quadrature2x2::Nq values in the output
 * vectors `vals` and `Dv`. Each entry of `Dv` is an array of three
 * numbers:
 * @f[ \left[
 * \frac{du}{dx}, \frac{dv}{dy}, \frac{1}{2}\left(\frac{du}{dy}+\frac{dv}{dx}\right)
 * \right] @f].
 */
void Quadrature_Vector::quadrature_point_values(const Vector2 *x, Vector2 *vals, double (*Dv)[3]) {
  for (unsigned int q = 0; q < Nq; q++) {
    vals[q].u = 0;
    vals[q].v = 0;
    double *Dvq = Dv[q];
    Dvq[0] = 0;
    Dvq[1] = 0;
    Dvq[2] = 0;
    const Germ *test = m_germs[q];
    for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
      vals[q] += test[k].val * x[k];
      Dvq[0]  += test[k].dx * x[k].u;
      Dvq[1]  += test[k].dy * x[k].v;
      Dvq[2]  += 0.5*(test[k].dy * x[k].u + test[k].dx * x[k].v);
    }
  }
}

/*! @brief Compute the values and symmetric gradient at the quadrature points of a vector-valued
  finite-element function with element-local degrees of freedom `x_nodal`.*/
/*! There should be room for Quadrature2x2::Nq values in the output vectors `vals`, `dx`, and `dy`.
  Each element of `dx` is the derivative of the vector-valued finite-element function in the x direction,
  and similarly for `dy`.
*/
void Quadrature_Vector::quadrature_point_values(const Vector2 *x, Vector2 *vals, Vector2 *dx, Vector2 *dy) {
  for (unsigned int q = 0; q < Nq; q++) {
    vals[q].u = 0;
    vals[q].v = 0;
    dx[q].u   = 0;
    dx[q].v   = 0;
    dy[q].u   = 0;
    dy[q].v   = 0;
    const Germ *test = m_germs[q];
    for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
      vals[q] += test[k].val * x[k];
      dx[q]   += test[k].dx * x[k];
      dy[q]   += test[k].dy * x[k];
    }
  }
}

//! The quadrature points on the reference square @f$ x,y=\pm 1/\sqrt{3} @f$.
const double Quadrature2x2::quadPoints[Quadrature2x2::Nq][2] =
  {{-0.57735026918962573, -0.57735026918962573},
   { 0.57735026918962573, -0.57735026918962573},
   { 0.57735026918962573,  0.57735026918962573},
   {-0.57735026918962573,  0.57735026918962573}};

//! The weights w_i for gaussian quadrature on the reference element with these quadrature points
const double Quadrature2x2::quadWeights[Quadrature2x2::Nq]  = {1.0, 1.0, 1.0, 1.0};

DirichletData::DirichletData()
  : m_indices(NULL), m_weight(1.0) {
  for (unsigned int k = 0; k < ShapeQ1::Nk; ++k) {
    m_indices_e[k] = 0;
  }
}

DirichletData::~DirichletData() {
  finish(NULL);
}

void DirichletData::init(const IceModelVec2Int *indices,
                         const IceModelVec *values,
                         double weight) {
  m_weight  = weight;

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
    MPI_Comm com = m_indices->get_grid()->ctx()->com();
    try {
      m_indices->end_access();
      m_indices = NULL;
    } catch (...) {
      handle_fatal_errors(com);
    }
  }

  if (values != NULL) {
    MPI_Comm com = values->get_grid()->ctx()->com();
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
  for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
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
  for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i = 0, j = 0;
      element.local_to_global(k, i, j);
      x_nodal[k] = (*m_values)(i, j);
    }
  }
}

void DirichletData_Scalar::enforce_homogeneous(const ElementMap &element, double* x_nodal) {
  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_nodal[k] = 0.;
    }
  }
}

void DirichletData_Scalar::fix_residual(double const *const *const x_global, double **r_global) {
  assert(m_values != NULL);

  const IceGrid &grid = *m_indices->get_grid();

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
  const IceGrid &grid = *m_indices->get_grid();

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
  const IceGrid &grid = *m_indices->get_grid();

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
  for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i = 0, j = 0;
      element.local_to_global(k, i, j);
      x_nodal[k] = (*m_values)(i, j);
    }
  }
}

void DirichletData_Vector::enforce_homogeneous(const ElementMap &element, Vector2* x_nodal) {
  element.nodal_values(*m_indices, m_indices_e);
  for (unsigned int k = 0; k < ShapeQ1::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_nodal[k].u = 0.0;
      x_nodal[k].v = 0.0;
    }
  }
}

void DirichletData_Vector::fix_residual(Vector2 const *const *const x_global, Vector2 **r_global) {
  assert(m_values != NULL);

  const IceGrid &grid = *m_indices->get_grid();

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
  const IceGrid &grid = *m_indices->get_grid();

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
  const IceGrid &grid = *m_indices->get_grid();

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

BoundaryQuadrature2::BoundaryQuadrature2(double dx, double dy) {

  const double jacobian_x = 0.5*dx;
  const double jacobian_y = 0.5*dy;

  // Note that all quadrature weights are 1.0 (and so they are implicitly included below).
  //
  // bottom
  m_weighted_jacobian[0] = jacobian_x;
  // right
  m_weighted_jacobian[1] = jacobian_y;
  // top
  m_weighted_jacobian[2] = jacobian_x;
  // left
  m_weighted_jacobian[3] = jacobian_y;

  const double C = 1.0 / sqrt(3);
  const double pts[n_sides][Nq][2] = {
    {{  -C, -1.0}, {   C, -1.0}}, // South
    {{ 1.0,   -C}, { 1.0,    C}}, // East
    {{  -C,  1.0}, {   C,  1.0}}, // North
    {{-1.0,   -C}, {-1.0,    C}}  // West
  };

  memset(m_germs, 0, n_sides*Nq*ShapeQ1::Nk*sizeof(Germ));

  for (unsigned int side = 0; side < n_sides; ++side) {
    for (unsigned int q = 0; q < Nq; ++q) {
      const double xi = pts[side][q][0];
      const double eta = pts[side][q][1];
      for (unsigned int k = 0; k < ShapeQ1::Nk; ++k) {
        m_germs[side][q][k] = ShapeQ1::eval(k, xi, eta);
        // convert from derivatives with respect to xi and eta to derivatives with respect to x and
        // y
        m_germs[side][q][k].dx /= jacobian_x;
        m_germs[side][q][k].dy /= jacobian_y;
      }
    }
  }
}

} // end of namespace fem
} // end of namespace pism
