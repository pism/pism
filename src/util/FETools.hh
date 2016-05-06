// Copyright (C) 2009--2011, 2013, 2014, 2015, 2016, 2017 Jed Brown nd Ed Bueler and Constantine Khroulev and David Maxwell
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

// The following are macros (instead of inline functions) so that error handling
// is less cluttered.  They should be replaced with empty macros when in
// optimized mode.

#ifndef _FETOOLS_H_
#define _FETOOLS_H_

#include <cassert>

#include <vector>

#include <petscmat.h>

#include "pism/util/Vector2.hh"

namespace pism {
class IceModelVec;
class IceModelVec2S;
class IceModelVec2Int;
class IceModelVec2V;
class IceGrid;
namespace fem {

//! Hard-wired maximum number of points a quadrature can use. This is used as the size of arrays
//! storing quadrature point values. Some of the entries in such an array may not be used.
static const unsigned int MAX_QUADRATURE_SIZE = 9;

//! \brief Classes for implementing the Finite Element Method on an IceGrid.
/*! \file FETools.hh We assume that the reader has a basic understanding of the finite element method. The
  following is a reminder of the method that also gives the background for the how to implement it
  on an IceGrid with the tools in this module.

  The IceGrid domain \f$\Omega\f$ is decomposed into a grid of rectangular physical elements indexed
  by indices (i,j):

  ~~~
  (0,1)       (1,1)
      ---------
      |   |   |
      ---------
      |   |   |
      ---------
  (0,0)       (1,0)
  ~~~

  The index of an element corresponds with the index of its lower-left vertex in the grid.

  The reference element is the square \f$[-1,1]\times[-1,1]\f$. For each physical element
  \f$E_{ij}\f$, there is an map \f$F_{ij}\f$ from the reference element \f$R\f$ to \f$E_{ij}\f$. In
  this implementation, the rectangles in the domain are all congruent, and the maps F_{ij} are all
  the same up to a translation.

  On the reference element, vertices are ordered as follows:

  ~~~
  3 o---------o 2
    |         |
    |         |
    |         |
  0 o---------o 1
  ~~~

  For each vertex \f$k\in\{0,1,2,3\}\f$ there is an element basis function \f$\phi_k\f$ that is bilinear, equals 1 at
  vertex \f$k\f$, and equals 0 at the remaining vertices.

  For each node \f$(i,j)\f$ in the physical domain there is a basis function \f$\psi_{ij}\f$ that equals 1 at
  vertex \f$(i,j)\f$, and equals zero at all other vertices, and that on element \f$E_{i'j'}\f$
  can be written as \f$\phi_k\circ F_{i',j'}^{-1}\f$ for some index \f$k\f$:
  \f[
  \psi_{ij}\Big|_{E_{i'j'}} = \phi_k\circ F_{i',j'}^{-1}.
  \f]

  The hat functions \f$\psi_{ij}\f$ form a basis of the function space of piecewise-linear functions.
  A (scalar) piecewise-linear function \f$v=v(x,y)\f$ on the domain is a linear combination
  \f[
  v = \sum_{i,j} c_{ij}\psi_{ij}.
  \f]

  Let \f$G(w,v)\f$ denote the weak form of the PDE under consideration. For example, for the
  scalar Poisson equation \f$-\Delta w = f\f$,
  \f[
  G(w,v) = \int_\Omega \nabla w \cdot \nabla v -f v \;dx.
  \f]
  The SSA weak form is more complicated, in particular because there is dimension 2 at each vertex,
  corresponding to x- and y-velocity components, but it is in many ways like the Poisson weak form.

  In the continuous problem we seek to find a trial function \f$w\f$ such that \f$G(w,v)=0\f$ for
  all suitable test functions \f$v\f$.

  In the discrete problem, we seek a finite element function
  \f$w_h\f$ such that \f$G(w_h,\psi_{ij})=0\f$ for all suitable indices \f$(i,j)\f$. For realistic
  problems, the integral defining \f$G\f$ cannot be evaluated exactly, but is approximated with some
  \f$G_h\f$ that arises from numerical quadrature rule: integration on an element \f$E\f$ is
  approximated with
  \f[
  \int_E f dx \approx \sum_{q} f(x_q) j_q
  \f]
  for certain points \f$x_q\f$ and weights \f$j_q\f$ (specific details are found in Quadrature).

  The unknown \f$w_h\f$ is represented by an IceVec, \f$w_h=\sum_{ij} c_{ij} \psi_{ij}\f$ where
  \f$c_{ij}\f$ are the coefficients of the vector. The solution of the finite element problem
  requires the following computations:

  -# Evaluation of the residuals \f$r_{ij} = G_h(w_h,\psi_{ij})\f$
  -# Evaluation of the Jacobian matrix
  \f[
  J_{(ij)\;(kl)}=\frac{d}{dc_{kl}}  G_h(w_h,\psi_{ij}).
  \f]

  Computations of these 'integrals' are done by adding up the contributions from each element (an
  ElementIterator helps with iterating over the elements). For a fixed element, there is a locally
  defined trial function \f$\hat w_h\f$ (with 4 degrees of freedom in the scalar case) and 4 local
  test functions \f$\hat\psi_k\f$, one for each vertex.

  The contribution to the integrals proceeds as follows (for concreteness
  in the case of computing the residual):

  - Extract from the global degrees of freedom \f$c\f$ defining \f$w_h\f$ the local degrees of
    freedom \f$d\f$ defining \f$\hat w_h\f$. (ElementMap)

  - Evaluate the local trial function \f$w_h\f$ (values and derivatives as needed) at the quadrature
    points \f$x_q\f$ (Quadrature)

  - Evaluate the local test functions (again values and derivatives) at the quadrature points.
    (Quadrature)

  - Obtain the quadrature weights \f$j_q\f$ for the element (Quadrature)

  - Compute the values of the integrand \f$G(\hat w_h,\psi_k)\f$ at each quadrature point (call
    these \f$g_{qk}\f$) and form the weighted sums \f$y_k=\sum_{q} j_q g_{qk}\f$.

  - Each sum \f$y_k\f$ is the contribution of the current element to a residual entry \f$r_{ij}\f$,
    where local degree of freedom \f$k\f$ corresponds with global degree of freedom \f$(i,j)\f$. The
    local contibutions now need to be added to the global residual vector (ElementMap).

  Computation of the Jacobian is similar, except that there are now multiple integrals per element
  (one for each local degree of freedom of \f$\hat w_h\f$).

  All of the above is a simplified description of what happens in practice. The complications below
  treated by the following classes, and discussed in their documentation:

  - Ghost elements (as well as periodic boundary conditions): ElementIterator
  - Dirichlet data: ElementMap
  - Vector valued functions: (ElementMap, Quadrature)

  The classes in this module are not intended to be a fancy finite element package. Their purpose is
  to clarify the steps that occur in the computation of residuals and Jacobians in SSAFEM, and to
  isolate and consolidate the hard steps so that they are not scattered about the code.
*/

//! Struct for gathering the value and derivative of a function at a point.
/*! Germ in meant in the mathematical sense, sort of. */
struct Germ {
  //! Function value.
  double val;
  //! Function deriviative with respect to x.
  double dx;
  //! Function derivative with respect to y.
  double dy;
};

//! Coordinates of a quadrature point, in the (xi, eta) coordinate space (i.e. on the reference
//! element).
struct QuadPoint {
  double xi;
  double eta;
};

//! Q0 element information.
// FIXME: not sure if Q0 is the right notation here.
namespace q0 {
//! Number of shape functions on a Q0 element.
const int n_chi = 4;
//! Evaluate a piecewise-constant shape function and its derivatives.
Germ chi(unsigned int k, const QuadPoint &p);
} // end of namespace q0

class ElementGeometry {
public:
  unsigned int n_sides() const;
  unsigned int incident_node(unsigned int side, unsigned int k) const;
  Vector2 normal(unsigned int side) const;
protected:
  ElementGeometry(unsigned int sides);

  virtual unsigned int incident_node_impl(unsigned int side, unsigned int k) const = 0;
  std::vector<Vector2> m_normals;
private:
  unsigned int m_n_sides;
};

class BoundaryQuadrature {
public:
  BoundaryQuadrature(unsigned int size);
  virtual ~BoundaryQuadrature();

  unsigned int n() const;

  double weight(unsigned int side, unsigned int q) const;

  const Germ& germ(unsigned int side, unsigned int q, unsigned int test_function) const;
protected:
  virtual double weight_impl(unsigned int side,
                             unsigned int q) const = 0;
  virtual const Germ& germ_impl(unsigned int side,
                                unsigned int q,
                                unsigned int test_function) const = 0;
private:
  unsigned int m_Nq;
};

//! Q1 element information.
namespace q1 {

//! Number of shape functions on a Q1 element.
const int n_chi = 4;
//! Number of sides per element.
const int n_sides = 4;
//! Evaluate a Q1 shape function and its derivatives with respect to xi and eta.
Germ chi(unsigned int k, const QuadPoint &p);

//! Nodes incident to a side. Used to extract nodal values and add contributions.

class Q1ElementGeometry : public ElementGeometry {
public:
  Q1ElementGeometry();
private:
  unsigned int incident_node_impl(unsigned int side, unsigned int k) const;
};


//! 2-point quadrature for sides of Q1 quadrilateral elements.
class BoundaryQuadrature2 : public BoundaryQuadrature {
public:
  BoundaryQuadrature2(double dx, double dy, double L);

protected:
  double weight_impl(unsigned int side, unsigned int q) const;
  const Germ& germ_impl(unsigned int side, unsigned int q, unsigned int test_function) const;
private:
  //! Number of quadrature points per side.
  static const unsigned int m_size = 2;
  double m_W[n_sides][m_size];
  Germ m_germs[n_sides][m_size][q1::n_chi];
};

} // end of namespace q1

//! @brief P1 element embedded in a Q1 element.
//! function at node 2).
namespace p1 {
//! Evaluate a P1 shape function and its derivatives with respect to xi and eta.
Germ chi(unsigned int k, const QuadPoint &p);

//! Number of sides per element.
const int n_sides = 3;

class P1ElementGeometry : public ElementGeometry {
public:
  P1ElementGeometry(unsigned int type, double dx, double dy);
private:
  unsigned int incident_node_impl(unsigned int side, unsigned int k) const;
  unsigned int m_type;
};

} // end of namespace p1


// define Germs, which is an array of q1::n_chi "Germ"s
typedef Germ Germs[q1::n_chi];

//! The mapping from global to local degrees of freedom.
/*! Computations of residual and Jacobian entries in the finite element method are done by iterating
  of the elements and adding the various contributions from each element. To do this, degrees of
  freedom from the global vector of unknowns must be remapped to element-local degrees of freedom
  for the purposes of local computation, (and the results added back again to the global residual
  and Jacobian arrays).

  An ElementMap mediates the transfer between element-local and global degrees of freedom. In this
  very concrete implementation, the global degrees of freedom are either scalars (double's) or
  vectors (Vector2's), one per node in the IceGrid, and the local degrees of freedom on the element
  are q1::n_chi (%i.e. four) scalars or vectors, one for each vertex of the element.

  The ElementMap is also (perhaps awkwardly) overloaded to mediate transfering locally computed
  contributions to residual and Jacobian matrices to their global counterparts.

  See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
*/
class ElementMap {
public:
  ElementMap(const IceGrid &grid);
  ~ElementMap();

  /*! @brief Extract nodal values for the element (`i`,`j`) from global array `x_global`
      into the element-local array `result`.
  */
  template<typename T>
  void nodal_values(T const* const* x_global, T* result) const {
    for (unsigned int k = 0; k < q1::n_chi; ++k) {
      int i = 0, j = 0;
      local_to_global(k, i, j);
      result[k] = x_global[j][i];   // note the indexing order
    }
  }

  /*! @brief Extract nodal values for the element (`i`,`j`) from global IceModelVec `x_global`
      into the element-local array `result`.
  */
  template<class C, typename T>
  void nodal_values(const C& x_global, T* result) const {
    for (unsigned int k = 0; k < q1::n_chi; ++k) {
      int i = 0, j = 0;
      local_to_global(k, i, j);
      result[k] = x_global(i, j);   // note the indexing order
    }
  }

  /*! @brief Get nodal values of an integer mask. */
  void nodal_values(const IceModelVec2Int &x_global, int *result) const;

  /*! @brief Add the values of element-local contributions `y` to the global vector `y_global`. */
  /*! The element-local array `local` should be an array of Nk values.
   *
   * Use this to add residual contributions.
   */
  template<typename T>
  void add_contribution(const T *local, T** y_global) const {
    for (unsigned int k = 0; k < fem::q1::n_chi; k++) {
      if (m_row[k].k == 1) {
        // skip rows marked as "invalid"
        continue;
      }
      const int
        i = m_row[k].i,
        j = m_row[k].j;
      y_global[j][i] += local[k];   // note the indexing order
    }
  }

  template<class C, typename T>
  void add_contribution(const T *local, C& y_global) const {
    for (unsigned int k = 0; k < fem::q1::n_chi; k++) {
      if (m_row[k].k == 1) {
        // skip rows marked as "invalid"
        continue;
      }
      const int
        i = m_row[k].i,
        j = m_row[k].j;
      y_global(i, j) += local[k];   // note the indexing order
    }
  }

  void reset(int i, int j);

  void mark_row_invalid(int k);
  void mark_col_invalid(int k);

  //! Convert a local degree of freedom index `k` to a global degree of freedom index (`i`,`j`).
  void local_to_global(int k, int &i, int &j) const {
    i = m_i + m_i_offset[k];
    j = m_j + m_j_offset[k];
  }

  /*! @brief Add Jacobian contributions. */
  void add_contribution(const double *K, Mat J) const;

private:
  //! Constant for marking invalid row/columns.
  //!
  //! Has to be negative because MatSetValuesBlockedStencil supposedly ignores negative indexes.
  //! Seems like it has to be negative and below the smallest allowed index, i.e. -2 and below with
  //! the stencil of width 1 (-1 *is* an allowed index). We use -2^30 and *don't* use PETSC_MIN_INT,
  //! because PETSC_MIN_INT depends on PETSc's configuration flags.
  static const int m_invalid_dof = -1073741824;
  static const int m_i_offset[q1::n_chi];
  static const int m_j_offset[q1::n_chi];

  //! Indices of the current element.
  int m_i, m_j;

  //! Stencils for updating entries of the Jacobian matrix.
  MatStencil m_row[q1::n_chi], m_col[q1::n_chi];

  // the grid (for marking ghost DOFs as "invalid")
  const IceGrid &m_grid;
};

//! Manages iterating over element indices.
/*! When computing residuals and Jacobians, there is a loop over all the elements
  in the IceGrid, and computations are done on each element.  The IceGrid
  has an underlying Petsc DA, and our processor does not own all of the nodes in the grid.
  Therefore we should not perform computation on all of the elements.  In general,
  an element will have ghost (o) and real (*) vertices:
  \verbatim
  o---*---*---*---o
  |   |   |   |   |
  o---*---*---*---o
  |   |   |   |   |
  o---o---o---o---o
  \endverbatim
  The strategy is to do computations on this processor on every element that has
  a vertex that is owned by this processor.  But we only update entries in the
  global residual and Jacobian matrices if the corresponding row corresponds to a
  vertex that we own.  In the worst case, where each vertex of an element is owned by
  a different processor, the computations for that element will be repeated four times,
  once for each processor.

  This same strategy also correctly deals with periodic boundary conditions. The way Petsc deals
  with periodic boundaries can be thought of as using a kind of ghost. So the rule still works:
  compute on all elements containg a real vertex, but only update rows corresponding to that real
  vertex.

  The calculation of what elements to index over needs to account for ghosts and the
  presence or absense of periodic boundary conditions in the IceGrid.  The ElementIterator performs
  that computation for you (see ElementIterator::xs and friends).

  See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
*/
class ElementIterator {
public:
  ElementIterator(const IceGrid &g);

  /*!\brief The total number of elements to be iterated over.  Useful for creating per-element storage.*/
  int element_count() {
    return xm*ym;
  }

  /*!\brief Convert an element index (`i`,`j`) into a flattened (1-d) array index, with the first
    element (`i`, `j`) to be iterated over corresponding to flattened index 0. */
  int flatten(int i, int j) {
    return (i-xs) + (j-ys)*xm;
  }

  //! x-coordinate of the first element to loop over.
  int xs;
  //! total number of elements to loop over in the x-direction.
  int xm;

  //! y-coordinate of the first element to loop over.
  int ys;
  //! total number of elements to loop over in the y-direction.
  int ym;

  //! x-index of the first local element.
  int lxs;
  //! total number local elements in x direction.
  int lxm;

  //! y-index of the first local element.
  int lys;
  //! total number local elements in y direction.
  int lym;
};

class Quadrature {
public:
  Quadrature(unsigned int N);
  virtual ~Quadrature();

  //! Quadrature size (the number of points).
  unsigned int n() const {
    return m_Nq;
  }

  //! Vaues of `q1::n_chi` trial functions and their derivatives with respect to `x` and `y`, for
  //! each of `n()` quadrature points.
  const Germs* test_function_values() const {
    return m_germs;
  }

  //! Determinant of the Jacobian of the map from the reference element to the physical element,
  //! evaluated at quadrature points and multiplied by corresponding quadrature weights.
  const double* weights() const {
    return m_W;
  }

  Germ test_function_values(unsigned int q, unsigned int k) const;
  double weights(unsigned int q) const;

protected:
  //! Number of quadrature points.
  const unsigned int m_Nq;

  double *m_W;

  Germs* m_germs;

  //! Jacobian of the map from the reference element to a physical element.
  double m_J[2][2];

  // pointer to a 2D shape function
  typedef Germ (*ShapeFunction2)(unsigned int k, const QuadPoint &p);

  void initialize(ShapeFunction2 f,
                  unsigned int n_chi,
                  const QuadPoint *points,
                  const double *weights);
};

//! Quadratures on a Q1 element with sides aligned along coordinate axes.
/**
  This code isolates the computation of the Jacobian of the map from the reference element to the
  physical element and its inverse, which are used to convert partial derivatives with respect to
  xi, eta to partial derivatives with respect to x and y.
 */
class UniformQxQuadrature : public Quadrature {
protected:
  UniformQxQuadrature(unsigned int size, double dx, double dy, double L);
};

//! Numerical integration of finite element functions.
/*! The core of the finite element method is the computation of integrals over elements.
  For nonlinear problems, or problems with non-constant coefficients (%i.e. any real problem)
  the integration has to be done approximately:
  \f[
  \int_E f(x)\; dx \approx \sum_q f(x_q) w_q
  \f]
  for certain quadrature points \f$x_q\f$ and weights \f$w_q\f$.  A quadrature is used
  to evaluate finite element functions at quadrature points, and to compute weights \f$w_q\f$
  for a given element.

  In this concrete implementation, the reference element \f$R\f$ is the square
  \f$[-1,1]\times[-1,1]\f$.  On a given element, nodes (o) and quadrature points (*)
  are ordered as follows:

  ~~~
  3 o------------------o  2
    |  3             2 |
    |    *        *    |
    |                  |
    |                  |
    |    *        *    |
    |  0            1  |
  0 o------------------o  1
  ~~~

  There are four quad points per element, which occur at \f$x,y=\pm 1/\sqrt{3}\f$. This corresponds
  to the tensor product of Gaussian integration on an interval that is exact for cubic functions on
  the interval.

  Integration on a physical element can be thought of as being done by change of variables. The
  quadrature weights need to be modified, and the Quadrature takes care of this. Because all
  elements in an IceGrid are congruent, the quadrature weights are the same for each element, and
  are computed upon initialization.

  See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
*/
class Q1Quadrature4 : public UniformQxQuadrature {
public:
  Q1Quadrature4(double dx, double dy, double L=1.0);
private:
  static const unsigned int m_size = 4;
};

//! The 9-point 2D Gaussian quadrature on the square [-1,1]*[-1,1].
class Q1Quadrature9 : public UniformQxQuadrature {
public:
  Q1Quadrature9(double dx, double dy, double L=1.0);
private:
  static const unsigned int m_size = 9;
};

//! The 16-point 2D Gaussian quadrature on the square [-1,1]*[-1,1].
class Q1Quadrature16 : public UniformQxQuadrature {
public:
  Q1Quadrature16(double dx, double dy, double L=1.0);
private:
  static const unsigned int m_size = 16;
};

class Q0Quadrature1e4 : public UniformQxQuadrature {
public:
  Q0Quadrature1e4(double dx, double dy, double L=1.0);
private:
  static const unsigned int m_size_1d = 100;
  static const unsigned int m_size = m_size_1d * m_size_1d;
};

class Q1Quadrature1e4 : public UniformQxQuadrature {
public:
  Q1Quadrature1e4(double dx, double dy, double L=1.0);
private:
  static const unsigned int m_size_1d = 100;
  static const unsigned int m_size = m_size_1d * m_size_1d;
};

//! Quadratures on a P1 element.
class P1Quadrature : public Quadrature {
protected:
  P1Quadrature(unsigned int size, unsigned int N,
               double dx, double dy, double L);
};

class P1Quadrature3 : public P1Quadrature {
public:
  P1Quadrature3(unsigned int N, double dx, double dy, double L);
private:
  static const unsigned int m_size = 3;
};

/*! @brief Compute the values at the quadrature points of a finite-element function.*/
/*! There should be room for Q.N() values in the output array `result`. */
template <typename T>
void quadrature_point_values(Quadrature &Q, const T *x, T *result) {
  const Germs *test = Q.test_function_values();
  const unsigned int n = Q.n();
  for (unsigned int q = 0; q < n; q++) {
    result[q] = 0.0;
    for (unsigned int k = 0; k < q1::n_chi; k++) {
      result[q] += test[q][k].val * x[k];
    }
  }
}

/*! @brief Compute the values and partial derivatives at the quadrature points of a finite-element
  function.*/
/*! There should be room for Q.N() values in the output vectors `vals`, `dx`, and `dy`.
  Each element of `dx` is the derivative of the vector-valued finite-element function in the x
  direction, and similarly for `dy`.
*/
template <typename T>
void quadrature_point_values(Quadrature &Q, const T *x, T *vals, T *dx, T *dy) {
  const Germs *test = Q.test_function_values();
  const unsigned int n = Q.n();
  for (unsigned int q = 0; q < n; q++) {
    vals[q] = 0.0;
    dx[q]   = 0.0;
    dy[q]   = 0.0;
    for (unsigned int k = 0; k < q1::n_chi; k++) {
      const Germ &psi = test[q][k];
      vals[q] += psi.val * x[k];
      dx[q]   += psi.dx  * x[k];
      dy[q]   += psi.dy  * x[k];
    }
  }
}

//* Parts shared by scalar and 2D vector Dirichlet data classes.
class DirichletData {
public:
  void constrain(ElementMap &element);
  operator bool() {
    return m_indices != NULL;
  }
protected:
  DirichletData();
  ~DirichletData();

  void init(const IceModelVec2Int *indices, const IceModelVec *values, double weight = 1.0);
  void finish(const IceModelVec *values);

  const IceModelVec2Int *m_indices;
  double m_indices_e[q1::n_chi];
  double m_weight;
};

class DirichletData_Scalar : public DirichletData {
public:
  DirichletData_Scalar(const IceModelVec2Int *indices, const IceModelVec2S *values,
                       double weight = 1.0);
  ~DirichletData_Scalar();

  void enforce(const ElementMap &element, double* x_e);
  void enforce_homogeneous(const ElementMap &element, double* x_e);
  void fix_residual(double const *const *const x_global, double **r_global);
  void fix_residual_homogeneous(double **r_global);
  void fix_jacobian(Mat J);
protected:
  const IceModelVec2S *m_values;
};

class DirichletData_Vector : public DirichletData {
public:
  DirichletData_Vector(const IceModelVec2Int *indices, const IceModelVec2V *values,
                       double weight);
  ~DirichletData_Vector();

  void enforce(const ElementMap &element, Vector2* x_e);
  void enforce_homogeneous(const ElementMap &element, Vector2* x_e);
  void fix_residual(Vector2 const *const *const x_global, Vector2 **r_global);
  void fix_residual_homogeneous(Vector2 **r);
  void fix_jacobian(Mat J);
protected:
  const IceModelVec2V *m_values;
};

} // end of namespace fem
} // end of namespace pism

#endif/* _FETOOLS_H_*/
