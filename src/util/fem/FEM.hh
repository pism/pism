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
#ifndef PISM_FEM_H
#define PISM_FEM_H


//! \brief Classes for implementing the Finite Element Method on an IceGrid.
/*! @file FEM.hh We assume that the reader has a basic understanding of the finite element method. The
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
    freedom \f$d\f$ defining \f$\hat w_h\f$. (Element)

  - Evaluate the local trial function \f$w_h\f$ (values and derivatives as needed) at the quadrature
    points \f$x_q\f$ (Quadrature)

  - Evaluate the local test functions (again values and derivatives) at the quadrature points.
    (Quadrature)

  - Obtain the quadrature weights \f$j_q\f$ for the element (Quadrature)

  - Compute the values of the integrand \f$G(\hat w_h,\psi_k)\f$ at each quadrature point (call
    these \f$g_{qk}\f$) and form the weighted sums \f$y_k=\sum_{q} j_q g_{qk}\f$.

  - Each sum \f$y_k\f$ is the contribution of the current element to a residual entry \f$r_{ij}\f$,
    where local degree of freedom \f$k\f$ corresponds with global degree of freedom \f$(i,j)\f$. The
    local contibutions now need to be added to the global residual vector (Element).

  Computation of the Jacobian is similar, except that there are now multiple integrals per element
  (one for each local degree of freedom of \f$\hat w_h\f$).

  All of the above is a simplified description of what happens in practice. The complications below
  treated by the following classes, and discussed in their documentation:

  - Ghost elements (as well as periodic boundary conditions): ElementIterator
  - Dirichlet data: Element
  - Vector valued functions: (Element, Quadrature)

  The classes in this module are not intended to be a fancy finite element package. Their purpose is
  to clarify the steps that occur in the computation of residuals and Jacobians in SSAFEM, and to
  isolate and consolidate the hard steps so that they are not scattered about the code.
*/

namespace pism {
namespace fem {

//! Struct for gathering the value and derivative of a function at a point.
/*! Germ in meant in the mathematical sense, sort of. */
struct Germ {
  //! Function value.
  double val;
  //! Function deriviative with respect to x.
  double dx;
  //! Function derivative with respect to y.
  double dy;
  //! Function derivative with respect to z.
  double dz;
};

//! Coordinates of a quadrature point, in the (xi, eta) coordinate space (i.e. on the
//! reference element).
struct QuadPoint {
  double xi;
  double eta;
  double zeta;
};

//! Hard-wired maximum number of points a quadrature can use. This is used as the size of
//! arrays storing quadrature point values. Some of the entries in such an array may not
//! be used.
const unsigned int MAX_QUADRATURE_SIZE = 9;

//! 1D (linear) elements
namespace linear {
const int n_chi = 2;
Germ chi(unsigned int k, const QuadPoint &pt);
} // end of namespace linear

//! Q0 element information
// FIXME: not sure if Q0 is the right notation here.
namespace q0 {
const int n_chi = 4;
const int n_sides = 4;
Germ chi(unsigned int k, const QuadPoint &p);
} // end of namespace q0

//! Q1 element information
namespace q1 {
const int n_chi = 4;
const int n_sides = 4;
Germ chi(unsigned int k, const QuadPoint &p);
} // end of namespace q1

//! P1 element information
namespace p1 {
Germ chi(unsigned int k, const QuadPoint &p);
const int n_chi = 3;
const int n_sides = 3;
} // end of namespace p1

enum ElementType {ELEMENT_Q = -1,
                  ELEMENT_P0 = 0, ELEMENT_P1 = 1, ELEMENT_P2 = 2, ELEMENT_P3 = 3,
                  ELEMENT_EXTERIOR};

ElementType element_type(int node_type[q1::n_chi]);

//! Q1 element information.
namespace q13d {

//! Number of shape functions on a Q1 element.
const int n_chi = 8;
//! Number of sides per element.
const int n_faces = 6;
//! Evaluate a Q1 shape function and its derivatives with respect to xi and eta.
Germ chi(unsigned int k, const QuadPoint &p);

/*! Nodes incident to a side. Used to extract nodal values and add contributions.
 *
 * The order of faces is used in Q1Element3Face::reset()
 */
const unsigned int incident_nodes[n_faces][4] =
  {{3, 0, 4, 7},                // 0 - left,   xi   = -1
   {1, 2, 6, 5},                // 1 - right,  xi   = +1
   {0, 1, 5, 4},                // 2 - front,  eta  = -1
   {2, 3, 7, 6},                // 3 - back,   eta  = +1
   {0, 1, 2, 3},                // 4 - bottom, zeta = -1
   {4, 5, 6, 7}                 // 5 - top,    zeta = +1
};

enum ElementFace {FACE_LEFT   = 0,
                  FACE_RIGHT  = 1,
                  FACE_FRONT  = 2,
                  FACE_BACK   = 3,
                  FACE_BOTTOM = 4,
                  FACE_TOP    = 5};

} // end of namespace q13d

} // end of namespace fem
} // end of namespace pism

#include "DirichletData.hh"

#include "Element.hh"

#include "ElementIterator.hh"

#include "Quadrature.hh"

#endif /* PISM_FEM_H */
