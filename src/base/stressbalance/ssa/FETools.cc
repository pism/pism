// Copyright (C) 2009--2011, 2013, 2014 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#include "FETools.hh"
#include "flowlaws.hh"
#include "IceGrid.hh"
#include <assert.h>

namespace pism {

const FEShapeQ1::ShapeFunctionSpec FEShapeQ1::shapeFunction[FEShapeQ1::Nk] =
  {FEShapeQ1::shape0, FEShapeQ1::shape1, FEShapeQ1::shape2, FEShapeQ1::shape3};

FEElementMap::FEElementMap(const IceGrid &g) {
  // Start by assuming ghost elements exist in all directions.
  // Elements are indexed by their lower left vertex.  If there is a ghost
  // element on the right, its i-index will be the same as the maximum
  // i-index of a non-ghost vertex in the local grid.
  xs = g.xs - 1;                    // Start at ghost to the left.
  int xf = g.xs + g.xm - 1; // End at ghost to the right.
  ys = g.ys - 1;                    // Start at ghost at the bottom.
  int yf = g.ys + g.ym - 1; // End at ghost at the top.

  lxs = g.xs;
  int lxf = lxs + g.xm - 1;
  lys = g.ys;
  int lyf = lys + g.ym - 1;

  // Now correct if needed. The only way there will not be ghosts is if the
  // grid is not periodic and we are up against the grid boundary.

  if (!(g.periodicity & X_PERIODIC))
    {
      // Leftmost element has x-index 0.
      if (xs < 0) {
        xs = 0;
      }
      // Rightmost vertex has index g.Mx-1, so the rightmost element has index g.Mx-2
      if (xf > g.Mx - 2)
        {
          xf = g.Mx - 2;
          lxf = g.Mx - 2;
        }
    }

  if (!(g.periodicity & Y_PERIODIC))
    {
      // Bottom element has y-index 0.
      if (ys < 0) {
        ys = 0;
      }
      // Topmost vertex has index g.My - 1, so the topmost element has index g.My - 2
      if (yf > g.My - 2)
        {
          yf = g.My - 2;
          lyf = g.My - 2;
        }
    }

  // Tally up the number of elements in each direction
  xm = xf - xs + 1;
  ym = yf - ys + 1;
  lxm = lxf - lxs + 1;
  lym = lyf - lys + 1;

}

FEDOFMap::FEDOFMap() {
  m_i = 0;
  m_j = 0;
  PetscMemzero(m_row, Nk*sizeof(MatStencil));
  PetscMemzero(m_col, Nk*sizeof(MatStencil));
}

FEDOFMap::~FEDOFMap() {
  // empty
}


/*! @brief Extract local degrees of freedom for element (`i`,`j`) from global vector `x_global` to
  local vector `x_local` (scalar-valued DOF version). */
void FEDOFMap::extractLocalDOFs(int i, int j, double const*const*x_global, double *x_local) const
{
  x_local[0] = x_global[i][j];
  x_local[1] = x_global[i + 1][j];
  x_local[2] = x_global[i + 1][j + 1];
  x_local[3] = x_global[i][j + 1];
}

void FEDOFMap::extractLocalDOFs(int i, int j, IceModelVec2S &x_global, double *x_local) const
{
  x_local[0] = x_global(i, j);
  x_local[1] = x_global(i + 1, j);
  x_local[2] = x_global(i + 1, j + 1);
  x_local[3] = x_global(i, j + 1);
}

/*! @brief Extract local degrees of freedom for element (`i`,`j`) from global vector `x_global` to
  local vector `x_local` (vector-valued DOF version).
*/
void FEDOFMap::extractLocalDOFs(int i, int j, Vector2 const*const*x_global,
                                Vector2 *x_local) const
{
  x_local[0] = x_global[i][j];
  x_local[1] = x_global[i + 1][j];
  x_local[2] = x_global[i + 1][j + 1];
  x_local[3] = x_global[i][j + 1];
}

void FEDOFMap::extractLocalDOFs(int i, int j, IceModelVec2V &x_global,
                                Vector2 *x_local) const
{
  x_local[0] = x_global(i, j);
  x_local[1] = x_global(i + 1, j);
  x_local[2] = x_global(i + 1, j + 1);
  x_local[3] = x_global(i, j + 1);
}


//! Extract scalar degrees of freedom for the element specified previously with FEDOFMap::reset
void FEDOFMap::extractLocalDOFs(double const*const*x_global, double *x_local) const
{
  extractLocalDOFs(m_i, m_j, x_global, x_local);
}

void FEDOFMap::extractLocalDOFs(IceModelVec2S &x_global, double *x_local) const
{
  extractLocalDOFs(m_i, m_j, x_global, x_local);
}

void FEDOFMap::extractLocalDOFs(Vector2 const*const*x_global, Vector2 *x_local) const
{
  extractLocalDOFs(m_i, m_j, x_global, x_local);
}

//! Extract vector degrees of freedom for the element specified previously with FEDOFMap::reset
void FEDOFMap::extractLocalDOFs(IceModelVec2V &x_global, Vector2 *x_local) const
{
  extractLocalDOFs(m_i, m_j, x_global, x_local);
}

//! Convert a local degree of freedom index `k` to a global degree of freedom index (`i`,`j`).
void FEDOFMap::localToGlobal(int k, int *i, int *j) {
  *i = m_i + kIOffset[k];
  *j = m_j + kJOffset[k];
}

/*!@brief Initialize the FEDOFMap to element (`i`, `j`) for the purposes of inserting into
  global residual and Jacobian arrays. */
void FEDOFMap::reset(int i, int j, const IceGrid &grid) {
  m_i = i; m_j = j;
  // The meaning of i and j for a PISM IceGrid and for a Petsc DA are swapped (the so-called
  // fundamental transpose.  The interface between PISM and Petsc is the stencils, so all
  // interactions with the stencils involve a transpose.
  m_col[0].i = j;
  m_col[1].i = j;
  m_col[2].i = j + 1;
  m_col[3].i = j + 1;

  m_col[0].j = i;
  m_col[1].j = i + 1;
  m_col[2].j = i + 1;
  m_col[3].j = i;

  memcpy(m_row, m_col, Nk*sizeof(m_col[0]));

  // We do not ever sum into rows that are not owned by the local rank.
  for (int k = 0; k < Nk; k++) {
    int pism_i = m_row[k].j, pism_j = m_row[k].i;
    if (pism_i < grid.xs || grid.xs + grid.xm - 1 < pism_i ||
        pism_j < grid.ys || grid.ys + grid.ym - 1 < pism_j) {
      markRowInvalid(k);
    }
  }
}

/*!@brief Mark that the row corresponding to local degree of freedom `k` should not be updated
  when inserting into the global residual or Jacobian arrays. */
void FEDOFMap::markRowInvalid(int k) {
  m_row[k].i = m_row[k].j = kDofInvalid;
}

/*!@brief Mark that the column corresponding to local degree of freedom `k` should not be updated
  when inserting into the global Jacobian arrays. */
void FEDOFMap::markColInvalid(int k) {
  m_col[k].i = m_col[k].j = kDofInvalid;
}

/*!@brief Add the values of element-local residual contributions `y` to the global residual
  vector `yg`. */
/*! The element-local residual should be an array of Nk values.*/
void FEDOFMap::addLocalResidualBlock(const Vector2 *y, Vector2 **yg) {
  for (unsigned int k = 0; k < Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    yg[m_row[k].j][m_row[k].i].u += y[k].u;
    yg[m_row[k].j][m_row[k].i].v += y[k].v;
  }
}

void FEDOFMap::addLocalResidualBlock(const double *y, double **yg) {
  for (unsigned int k = 0; k < Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    yg[m_row[k].j][m_row[k].i] += y[k];
  }
}

void FEDOFMap::addLocalResidualBlock(const Vector2 *y, IceModelVec2V &y_global) {
  for (unsigned int k = 0; k < Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    y_global(m_row[k].j, m_row[k].i).u += y[k].u;
    y_global(m_row[k].j, m_row[k].i).v += y[k].v;
  }
}

void FEDOFMap::addLocalResidualBlock(const double *y, IceModelVec2S &y_global) {
  for (unsigned int k = 0; k < Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    y_global(m_row[k].j, m_row[k].i) += y[k];
  }
}

//! Add the contributions of an element-local Jacobian to the global Jacobian vector.
/*! The element-local Jacobian should be givnen as a row-major array of Nk*Nk values in the
  scalar case or (2Nk)*(2Nk) values in the vector valued case. */
PetscErrorCode FEDOFMap::addLocalJacobianBlock(const double *K, Mat J) {
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J, Nk, m_row, Nk, m_col, K, ADD_VALUES); CHKERRQ(ierr);
  return 0;
}

//! Set a diagonal entry for global degree of freedom (`i` ,`j`) in a Jacobian matrix
/*! This is an unhappy hack for supporting Dirichlet constrained degrees of freedom.
  In the scalar valued case, `K` should point to a single value, and in the vector case,
  it should point to 4 (=2x2) values for the (2x2) block correspoinding to to u-u, u-v, v-u, and v-v
  interactions at grid point (`i`, `j`).  Sheesh.*/
PetscErrorCode FEDOFMap::setJacobianDiag(int i, int j, const double*K, Mat J) {
  MatStencil row;
  row.i = j; row.j = i;
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, K, INSERT_VALUES); CHKERRQ(ierr);
  return 0;
}

const int FEDOFMap::kIOffset[4] = {0, 1, 1, 0};
const int FEDOFMap::kJOffset[4] = {0, 0, 1, 1};


FEQuadrature_Scalar::FEQuadrature_Scalar() {
  PetscMemzero(m_tmp, Nk*sizeof(double));
}

//! Obtain the weights @f$ w_q @f$ for quadrature.
void FEQuadrature::getWeightedJacobian(double *jxw) {
  for (int q=0; q<Nq; q++)
    {
      jxw[q] = m_jacobianDet * quadWeights[q];
    }
}

//! Obtain the weights @f$ w_q @f$ for quadrature.
void FEQuadrature::init(const IceGrid &grid, double L) {
  // Since we use uniform cartesian coordinates, the Jacobian is
  // constant and diagonal on every element.
  //
  // Note that the reference element is @f$ [-1,1]^2 @f$ hence the
  // extra factor of 1/2.
  double jacobian_x = 0.5*grid.dx / L;
  double jacobian_y = 0.5*grid.dy / L;
  m_jacobianDet = jacobian_x*jacobian_y;

  FEShapeQ1 shape;
  for (int q = 0; q < Nq; q++) {
    for (int k = 0; k < Nk; k++) {
      shape.eval(k, quadPoints[q][0], quadPoints[q][1], &m_germs[q][k]);
      m_germs[q][k].dx /= jacobian_x;
      m_germs[q][k].dy /= jacobian_y;
    }
  }
}

FEQuadrature_Vector::FEQuadrature_Vector() {
  PetscMemzero(m_tmp, Nk*sizeof(Vector2));
}

//! Return the values at all quadrature points of all shape functions.
//* The return value is an Nq by Nk array of FEFunctionGerms. */
const FEQuadrature::FEFunctionGermArray* FEQuadrature::testFunctionValues()
{
  return m_germs;
}

//! Return the values of all shape functions at quadrature point `q`
//* The return value is an array of Nk FEFunctionGerms. */
const FEFunctionGerm *FEQuadrature::testFunctionValues(int q) {
  return m_germs[q];
}

//! Return the values at quadrature point `q` of shape function `k`.
const FEFunctionGerm *FEQuadrature::testFunctionValues(int q, int k) {
  return m_germs[q] + k;
}


/*! @brief Compute the values at the quadrature ponits of a scalar-valued
  finite-element function with element-local degrees of freedom `x_local`.*/
/*! There should be room for FEQuadrature::Nq values in the output vector `vals`. */
void FEQuadrature_Scalar::computeTrialFunctionValues(const double *x_local, double *vals) {
  for (int q = 0; q < Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0;
    for (int k = 0; k < Nk; k++) {
      vals[q] += test[k].val * x_local[k];
    }
  }
}

/*! @brief Compute the values and first derivatives at the quadrature
  points of a scalar-valued finite-element function with element-local
  degrees of freedom `x_local`.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors `vals`, `dx`,
  and `dy`. */
void FEQuadrature_Scalar::computeTrialFunctionValues(const double *x_local, double *vals, double *dx, double *dy) {
  for (int q = 0; q < Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0; dx[q] = 0; dy[q] = 0;
    for (int k = 0; k < Nk; k++) {
      vals[q] += test[k].val * x_local[k];
      dx[q]   += test[k].dx * x_local[k];
      dy[q]   += test[k].dy * x_local[k];
    }
  }
}

/*! @brief Compute the values at the quadrature points on element (`i`,`j`)
  of a scalar-valued finite-element function with global degrees of freedom `x`.*/
/*! There should be room for FEQuadrature::Nq values in the output vector `vals`. */
void FEQuadrature_Scalar::computeTrialFunctionValues(int i, int j, const FEDOFMap &dof,
                                                     double const*const*x_global, double *vals) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals);
}


void FEQuadrature_Scalar::computeTrialFunctionValues(int i, int j, const FEDOFMap &dof,
                                                     IceModelVec2S &x_global, double *vals) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals);
}

/*! @brief Compute the values and first derivatives at the quadrature
  points on element (`i`,`j`) of a scalar-valued finite-element function
  with global degrees of freedom `x`.*/
/*! There should be room for FEQuadrature::Nq values in the output
  vectors `vals`, `dx`, and `dy`. */
void FEQuadrature_Scalar::computeTrialFunctionValues(int i, int j,
                                                     const FEDOFMap &dof, double const*const*x_global,
                                                     double *vals, double *dx, double *dy) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals, dx, dy);
}

void FEQuadrature_Scalar::computeTrialFunctionValues(int i, int j,
                                                     const FEDOFMap &dof, IceModelVec2S &x_global,
                                                     double *vals, double *dx, double *dy) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals, dx, dy);
}

/*! @brief Compute the values at the quadrature points of a vector-valued
  finite-element function with element-local degrees of freedom `x_local`.*/
/*! There should be room for FEQuadrature::Nq values in the output vector `vals`. */
void FEQuadrature_Vector::computeTrialFunctionValues(const Vector2 *x_local, Vector2 *result) {
  for (int q = 0; q < Nq; q++) {
    result[q].u = 0;
    result[q].v = 0;
    const FEFunctionGerm *test = m_germs[q];
    for (int k = 0; k < Nk; k++) {
      result[q].u += test[k].val * x_local[k].u;
      result[q].v += test[k].val * x_local[k].v;
    }
  }
}

/*! @brief Compute the values and symmetric gradient at the quadrature
 *         points of a vector-valued finite-element function with
 *         element-local degrees of freedom `x_local`.
 *
 * There should be room for FEQuadrature::Nq values in the output
 * vectors `vals` and `Dv`. Each entry of `Dv` is an array of three
 * numbers:
 * @f[ \left[
 * \frac{du}{dx}, \frac{dv}{dy}, \frac{1}{2}\left(\frac{du}{dy}+\frac{dv}{dx}\right)
 * \right] @f].
 */
void FEQuadrature_Vector::computeTrialFunctionValues(const Vector2 *x_local, Vector2 *vals, double (*Dv)[3]) {
  for (int q = 0; q < Nq; q++) {
    vals[q].u = 0; vals[q].v = 0;
    double *Dvq = Dv[q];
    Dvq[0] = 0; Dvq[1] = 0; Dvq[2] = 0;
    const FEFunctionGerm *test = m_germs[q];
    for (int k = 0; k < Nk; k++) {
      vals[q].u += test[k].val * x_local[k].u;
      vals[q].v += test[k].val * x_local[k].v;
      Dvq[0] += test[k].dx * x_local[k].u;
      Dvq[1] += test[k].dy * x_local[k].v;
      Dvq[2] += 0.5*(test[k].dy*x_local[k].u + test[k].dx*x_local[k].v);
    }
  }
}

/*! @brief Compute the values and symmetric gradient at the quadrature points of a vector-valued
  finite-element function with element-local degrees of freedom `x_local`.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors `vals`, \ dx, and `dy`.
  Each element of `dx` is the derivative of the vector-valued finite-element function in the x direction,
  and similarly for `dy`.
*/
void FEQuadrature_Vector::computeTrialFunctionValues(const Vector2 *x_local, Vector2 *vals, Vector2 *dx, Vector2 *dy) {
  for (int q = 0; q < Nq; q++) {
    vals[q].u = 0; vals[q].v = 0;
    dx[q].u = 0; dx[q].v = 0;
    dy[q].u = 0; dy[q].v = 0;
    const FEFunctionGerm *test = m_germs[q];
    for (int k = 0; k < Nk; k++) {
      vals[q].u += test[k].val * x_local[k].u;
      vals[q].v += test[k].val * x_local[k].v;
      dx[q].u += test[k].dx * x_local[k].u;
      dx[q].v += test[k].dx * x_local[k].v;
      dy[q].u += test[k].dy * x_local[k].u;
      dy[q].v += test[k].dy * x_local[k].v;
    }
  }
}


/*! @brief Compute the values at the quadrature points of a vector-valued
  finite-element function on element (`i`,`j`) with global degrees of freedom `x_global`.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors `vals`. */
void FEQuadrature_Vector::computeTrialFunctionValues(int i, int j, const FEDOFMap &dof,
                                                     Vector2 const*const*x_global, Vector2 *vals) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals);
}

void FEQuadrature_Vector::computeTrialFunctionValues(int i, int j, const FEDOFMap &dof,
                                                     IceModelVec2V &x_global, Vector2 *vals) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals);
}

/*! @brief Compute the values and symmetric gradient at the quadrature points of a vector-valued
  finite-element function on element (`i`,`j`) with global degrees of freedom `x_global`.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors `vals` and `Dv`.
  Each entry of `Dv` is an array of three numbers:
  @f[\left[\frac{du}{dx},\frac{dv}{dy},\frac{1}{2}\left(\frac{du}{dy}+\frac{dv}{dx}\right)\right]@f].
*/
void FEQuadrature_Vector::computeTrialFunctionValues(int i, int j, const FEDOFMap &dof,
                                                     Vector2 const*const* x_global,
                                                     Vector2 *vals, double (*Dv)[3]) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals, Dv);
}

void FEQuadrature_Vector::computeTrialFunctionValues(int i, int j, const FEDOFMap &dof,
                                                     IceModelVec2V &x_global,
                                                     Vector2 *vals, double (*Dv)[3]) {
  dof.extractLocalDOFs(i, j, x_global, m_tmp);
  computeTrialFunctionValues(m_tmp, vals, Dv);
}

//! The quadrature points on the reference square @f$ x,y=\pm 1/\sqrt{3} @f$.
const double FEQuadrature::quadPoints[FEQuadrature::Nq][2] =
  {{-0.57735026918962573, -0.57735026918962573},
   { 0.57735026918962573, -0.57735026918962573},
   { 0.57735026918962573,  0.57735026918962573},
   {-0.57735026918962573,  0.57735026918962573}};

//! The weights w_i for gaussian quadrature on the reference element with these quadrature points
const double FEQuadrature::quadWeights[FEQuadrature::Nq]  = {1.0, 1.0, 1.0, 1.0};

DirichletData::DirichletData()
  : m_indices(NULL), m_weight(1.0) {
  for (int k = 0; k < FEQuadrature::Nk; ++k) {
    m_indices_e[k] = 0;
  }
}

DirichletData::~DirichletData() {
  if (m_indices != NULL) {
    PetscErrorCode ierr;
    ierr = verbPrintf(1, m_indices->get_grid()->com,
                      "Warning: DirichletData destructing with IceModelVecs still accessed."
                      " Looks like DirichletData::finish() was not called.");
    CHKERRCONTINUE(ierr);
  }
}

PetscErrorCode DirichletData::init(IceModelVec2Int *indices, double weight) {
  PetscErrorCode ierr = init_impl(indices, NULL, weight); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode DirichletData::finish() {
  PetscErrorCode ierr = finish_impl(NULL); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode DirichletData::init_impl(IceModelVec2Int *indices, IceModelVec *values,
                                        double weight) {
  PetscErrorCode ierr;
  m_weight  = weight;

  if (indices != NULL) {
    ierr = indices->begin_access(); CHKERRQ(ierr);
    m_indices = indices;
  }

  if (values != NULL) {
    ierr = values->begin_access(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode DirichletData::finish_impl(IceModelVec *values) {
  PetscErrorCode ierr;
  if (m_indices != NULL) {
    ierr = m_indices->end_access(); CHKERRQ(ierr);
    m_indices = NULL;
  }

  if (values != NULL) {
    ierr = values->end_access(); CHKERRQ(ierr);
  }

  return 0;
}

void DirichletData::constrain(FEDOFMap &dofmap) {
  dofmap.extractLocalDOFs(*m_indices, m_indices_e);
  for (int k = 0; k < FEQuadrature::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      // Mark any kind of Dirichlet node as not to be touched
      dofmap.markRowInvalid(k);
      dofmap.markColInvalid(k);
    }
  }
}

// Scalar version

DirichletData_Scalar::DirichletData_Scalar()
  : m_values(NULL) {
}

PetscErrorCode DirichletData_Scalar::init(IceModelVec2Int *indices, IceModelVec2S *values,
                                          double weight) {
  m_values  = values;
  PetscErrorCode ierr = init_impl(indices, m_values, weight); CHKERRQ(ierr);

  return 0;
}

void DirichletData_Scalar::update(FEDOFMap &dofmap, double* x_local) {
  assert(m_values != NULL);

  dofmap.extractLocalDOFs(*m_indices, m_indices_e);
  for (int k = 0; k < FEQuadrature::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i, j;
      dofmap.localToGlobal(k, &i, &j);
      x_local[k] = (*m_values)(i,j);
    }
  }
}

void DirichletData_Scalar::update_homogeneous(FEDOFMap &dofmap, double* x_local) {
  dofmap.extractLocalDOFs(*m_indices, m_indices_e);
  for (int k = 0; k < FEQuadrature::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_local[k] = 0.;
    }
  }
}

void DirichletData_Scalar::fix_residual(double **x_global, double **r_global) {
  assert(m_values != NULL);

  IceGrid &grid = *m_indices->get_grid();

  // For each node that we own:
  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {
      if ((*m_indices)(i, j) > 0.5) {
        // Enforce explicit dirichlet data.
        r_global[i][j] = m_weight * (x_global[i][j] - (*m_values)(i,j));
      }
    }
  }
}

void DirichletData_Scalar::fix_residual_homogeneous(double **r_global) {
  IceGrid *grid = m_indices->get_grid();

  // For each node that we own:
  for (int i = grid->xs; i < grid->xs + grid->xm; i++) {
    for (int j = grid->ys; j < grid->ys + grid->ym; j++) {
      if ((*m_indices)(i, j) > 0.5) {
        // Enforce explicit dirichlet data.
        r_global[i][j] = 0.0;
      }
    }
  }
}

PetscErrorCode DirichletData_Scalar::fix_jacobian(Mat J) {
  PetscErrorCode ierr;
  IceGrid *grid = m_indices->get_grid();

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.

  const double identity = m_weight;
  for (int i = grid->xs; i < grid->xs + grid->xm; i++) {
    for (int j = grid->ys; j < grid->ys + grid->ym; j++) {
      if ((*m_indices)(i, j) > 0.5) {
        MatStencil row;
        // Transpose shows up here!
        row.j = i; row.i = j;
        ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, &identity,
                                          ADD_VALUES); CHKERRQ(ierr);
      }
    }
  }
  return 0;
}

PetscErrorCode DirichletData_Scalar::finish() {
  PetscErrorCode ierr = finish_impl(m_values); CHKERRQ(ierr);
  m_values = NULL;
  return 0;
}

// Vector version

DirichletData_Vector::DirichletData_Vector()
  : m_values(NULL) {
}

PetscErrorCode DirichletData_Vector::init(IceModelVec2Int *indices, IceModelVec2V *values,
                                          double weight) {
  m_values = values;
  PetscErrorCode ierr = init_impl(indices, m_values, weight); CHKERRQ(ierr);

  return 0;
}

void DirichletData_Vector::update(FEDOFMap &dofmap, Vector2* x_local) {
  assert(m_values != NULL);

  dofmap.extractLocalDOFs(*m_indices, m_indices_e);
  for (int k = 0; k < FEQuadrature::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i, j;
      dofmap.localToGlobal(k, &i, &j);
      x_local[k].u = (*m_values)(i, j).u;
      x_local[k].v = (*m_values)(i, j).v;
    }
  }
}

void DirichletData_Vector::update_homogeneous(FEDOFMap &dofmap, Vector2* x_local) {
  dofmap.extractLocalDOFs(*m_indices, m_indices_e);
  for (int k = 0; k < FEQuadrature::Nk; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_local[k].u = 0.0;
      x_local[k].v = 0.0;
    }
  }
}

void DirichletData_Vector::fix_residual(Vector2 **x_global, Vector2 **r_global) {
  assert(m_values != NULL);

  IceGrid *grid = m_indices->get_grid();

  // For each node that we own:
  for (int i = grid->xs; i < grid->xs + grid->xm; i++) {
    for (int j = grid->ys; j < grid->ys + grid->ym; j++) {
      if ((*m_indices)(i, j) > 0.5) {
        // Enforce explicit dirichlet data.
        r_global[i][j].u = m_weight * (x_global[i][j].u - (*m_values)(i, j).u);
        r_global[i][j].v = m_weight * (x_global[i][j].v - (*m_values)(i, j).v);
      }
    }
  }
}

void DirichletData_Vector::fix_residual_homogeneous(Vector2 **r_global) {
  IceGrid &grid = *m_indices->get_grid();

  // For each node that we own:
  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {
      if ((*m_indices)(i, j) > 0.5) {
        // Enforce explicit dirichlet data.
        r_global[i][j].u = 0.0;
        r_global[i][j].v = 0.0;
      }
    }
  }
}

PetscErrorCode DirichletData_Vector::fix_jacobian(Mat J) {
  PetscErrorCode ierr;
  IceGrid &grid = *m_indices->get_grid();

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.

  const double identity[4] = {m_weight, 0,
                              0, m_weight};
  for (int i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (int j=grid.ys; j<grid.ys+grid.ym; j++) {
      if ((*m_indices)(i, j) > 0.5) {
        MatStencil row;
        // Transpose shows up here!
        row.j = i; row.i = j;
        ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, identity,
                                          ADD_VALUES); CHKERRQ(ierr);
      }
    }
  }
  return 0;
}

PetscErrorCode DirichletData_Vector::finish() {
  PetscErrorCode ierr = finish_impl(m_values); CHKERRQ(ierr);
  m_values = NULL;
  return 0;
}

} // end of namespace pism
