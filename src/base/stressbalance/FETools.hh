// Copyright (C) 2009--2011 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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


#include <petscmat.h>
#include "iceModelVec.hh"       // to get PISMVector2

//! \file 
//! \brief Classes for implementing the Finite Element Method on an IceGrid.
/*! \file
We assume that the reader has a basic understanding of the finite element method
at the level of ?????.  The following is a reminder of the method that also
gives the background for the how to implement it on an IceGrid with the tools in this 
module.

The IceGrid domain \f$\Omega\f$ is decomposed into a grid of rectangular physical elements indexed by indices (i,j):

\verbatim
  (0,1)   (1,1)
    ---------
    |   |   |
    ---------
    |   |   |
    ---------
   (0,0) (1,0)   
\endverbatim

The index of an element corresponds with the index of its lower-left vertex in the grid.

The reference element is the square \f$[-1,1]\times[-1,1]\f$. For each physical element \f$E_{ij}\f$, there is an 
map \f$F_{ij}\f$ from the reference element \f$R\f$ to \f$E_{ij}\f$.  
In this implementation, the rectangles in the domain are all congruent, and the maps F_{ij} are all the same up 
to a translation.

On the reference element, vertices are ordered as follows:

\verbatim
         3 o---------o  2
           |         |
           |         |
           |         |
        0  o---------o  1
\endverbatim

For each vertex \f$k\f$ there is an element basis function \f$\phi_k\f$ that is bilinear, equals 1 at
vertex \f$k\f$, and equals 0 at the remaining vertices.  

For each node \f$(i',j')\f$ in the physical domain there is a basis function that equals 1 at
vertex \f$(i',j')\f$, equals zero at all other vertices, and that on element \f$(i,j)\f$
can be written as \f$\phi_k\circ F_{i,j}^{-1}\f$ for some index \f$k\f$.

A (scalar) finite element function \f$f\f$ on the domain is then a linear combination 
\f[
f_h = \sum_{i,j} c_{ij}\psi_{ij}.
\f]

Let \f$G(w,\psi)\f$ denote the weak form of the PDE under consideration.  For 
example, for the scalar Poisson equation \f$-\Delta w = f\f$,
\f[
G(w,\psi) = \int_\Omega \nabla w \cdot \nabla \psi -f\psi\;dx.
\f]
In the continuous problem we seek to 
find a trial function \f$w\f$ such that \f$G(w,\psi)=0\f$ for all suitable test 
functions \f$\psi\f$.  In the discrete problem, we seek a finite element function \f$w_h\f$
such that \f$G(w_h,\psi_{ij})=0\f$ for all suitable indices \f$(i,j)\f$.
For realistice problems, the integral given by \f$G\f$ cannot be evaluated exactly, 
but is approximated with some \f$G_h\f$ that arises from numerical quadrature quadrature rule:
integration on an element \f$E\f$ is approximated with 
\f[
\int_E f dx \approx \sum_{q} f(x_q) w_q
\f]
for certain points \f$x_q\f$ and weights \f$j_q\f$ (specific details are found in FEQuadrature).

The unknown \f$w_h\f$ is represented by an IceVec, \f$w_h=\sum_{ij} c_{ij} \psi_{ij}\f$ where
\f$c_{ij}\f$ are the coefficients of the vector.  The solution of the finite element problem
requires the following computations:

-# Evaluation of the residuals \f$r_{ij} = G_h(w_h,\psi_{ij})\f$
-# Evaluation of the Jacobian matrix
\f[
J_{(ij)\;(kl)}=\frac{d}{dc_{kl}}  G_h(w_h,\psi_{ij}).
\f]

Computations of these 'integrals' are done by adding up the contributions
from each element (an FEElementMap helps with iterating over the elements).  
For a fixed element, there is a locally defined trial
function \f$\hat w_h\f$ (with 4 degrees of freedom in the scalar case)  and
4 local test functions \f$\hat\psi_k\f$, one for each vertex.  

The contribution to the integrals proceeds as follows (for concreteness
in the case of computing the residual):

- Extract from the global degrees of freedom \f$c\f$ defining \f$w_h\f$
the local degrees of freedom \f$d\f$ defining \f$\hat w_h\f$. (FEDOFMap)
- Evaluate the local trial function \f$w_h\f$ (values and derivatives as needed)
at the quadrature points \f$x_q\f$ (FEQuadrature)
- Evaluate the local test functions (again values and derivatives)
at the quadrature points. (FEQuadrature)
- Obtain the quadrature weights $j_q$ for the element (FEQuadrature)
- Compute the values of the integrand \f$g(\hat w_h,\psi_k)\f$
  at each quadrature point (call these \f$g_{qk}\f$) and
  form the weighted sums \f$y_k=\sum_{q} j_q g_{qk}\f$.
- Each sum \f$y_k\f$ is the contribution of the current element to
a residual entry \f$r_{ij}\f$, where local degree of freedom \f$k\f$
corresponds with global degree of freedom \f$(i,j)\f$. The
local contibutions now need to be added to the global residual vector (FEDOFMap).

Computation of the Jacobian is similar, except that there are now
multiple integrals per element (one for each local degree of freedom of
\f$\hat w_h\f$).

All of the above is a simplified description of what happens in practice.
The complications below treated by the following classes, and discussed
in their documentation:

- Ghost elements (as well as periodic boundary conditions): FEElementMap
- Dirichlet data: FEDOFMap
- Vector valued functions: (FEDOFMap, FEQuadrature)

The classes in this module are not intended
to be a fancy finite element package.
Their purpose is to clarify the steps that occur in the computation of
residuals and Jacobians in SSAFEM, and to isolate and consolodate 
the hard steps so that they are not scattered about the code.
*/

//! Struct for gathering the value and derivative of a function at a point.
/*! Germ in meant in the mathematical sense, sort of. */
struct FEFunctionGerm
{
  PetscReal val,  //!< Function value.
             dx,  //!< Function deriviative with respect to x.
             dy;  //!< Function derivative with respect to y.
};

//! Struct for gathering the value and derivative of a vector valued function at a point.
/*! Germ in meant in the mathematical sense, sort of. */
struct FEVector2Germ
{
  PISMVector2  val,  //!< Function value.
                dx,  //!< Function deriviative with respect to x.
                dy;  //!< Function derivative with respect to y.
};


//! Computation of Q1 shape function values.
/*! The Q1 shape functions are bilinear functions.  On a rectangular
element, there are four (FEShapeQ1::Nk) basis functions, each equal
to 1 at one vertex and equal to zero at the remainder.

The FEShapeQ1 class consolodates the computation of the values and
derivatives of these basis functions. */
class FEShapeQ1 {
public:
  virtual ~FEShapeQ1() {}

  //! Compute values and derivatives of the shape function supported at node 0
  static void shape0(PetscReal x, PetscReal y, FEFunctionGerm *value)
  {
    value->val = (1-x)*(1-y)/4.;
    value->dx = -(1-y)/4.;
    value->dy = -(1-x)/4;
  }
  //! Compute values and derivatives of the shape function supported at node 1
  static void shape1(PetscReal x, PetscReal y, FEFunctionGerm *value)
  {
    value->val = (1+x)*(1-y)/4.;
    value->dx =  (1-y)/4.;
    value->dy = -(1+x)/4;
  }
  //! Compute values and derivatives of the shape function supported at node 2
  static void shape2(PetscReal x, PetscReal y, FEFunctionGerm *value)
  {
    value->val = (1+x)*(1+y)/4.;
    value->dx =  (1+y)/4.;
    value->dy =  (1+x)/4;
  }
  //! Compute values and derivatives of the shape function supported at node 3
  static void shape3(PetscReal x, PetscReal y, FEFunctionGerm *value)
  {
    value->val = (1-x)*(1+y)/4.;
    value->dx =  -(1+y)/4.;
    value->dy =   (1-x)/4;
  }

  //! The number of basis shape functions.
  static const PetscInt Nk = 4;

  //! A table of function pointers, one for each shape function.
  typedef void (*ShapeFunctionSpec)(PetscReal,PetscReal,FEFunctionGerm*);
  static const ShapeFunctionSpec shapeFunction[Nk];
  
  //! Evaluate shape function \a k at (\a x,\a y) with values returned in \a germ.
  virtual void eval(PetscInt k, PetscReal x, PetscReal y,FEFunctionGerm*germ){
    shapeFunction[k](x,y,germ);
  }
};


// Computes the closest integer to maskvalue, with integers of the form n/2 rounded up.
int PismIntMask(PetscScalar maskvalue);


//! The mapping from global to local degrees of freedom.
/*! Computations of residual and Jacobian entries in the finite element method are
done by iterating of the elements and adding the various contributions from each element.  
To do this, degrees of freedom from the global vector of unknowns must be remapped to
element-local degrees of freedom for the purposes of local computation, 
(and the results added back again to the global residual and Jacobian arrays).

An FEDOFMap mediates the transfer between element-local and global degrees of freedom.
In this very concrete implementation, the global degrees of freedom are either
scalars (PetscReal's) or vectors (PISMVector2's), one per node in the IceGrid, 
and the local degrees of freedom on the element are FEDOFMap::Nk (%i.e. four) scalars or vectors, one 
for each vertex of the element.

The FEDOFMap is also (perhaps awkwardly) overloaded to also mediate transfering locally
computed contributions to residual and Jacobian matricies to their global
counterparts.

See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
*/
class FEDOFMap
{
public:
  FEDOFMap() {};
  
  virtual ~FEDOFMap() {};

  virtual void extractLocalDOFs(PetscInt i, PetscInt j, PetscReal const*const*xg, PetscReal *x) const;
  virtual void extractLocalDOFs(PetscInt i, PetscInt j, PISMVector2 const*const*xg, PISMVector2 *x) const;

  virtual void extractLocalDOFs(PetscReal const*const*xg, PetscReal *x) const;
  virtual void extractLocalDOFs(PISMVector2 const*const*xg, PISMVector2 *x) const;

  virtual void reset(PetscInt i, PetscInt j, const IceGrid &g);
  
  virtual void markRowInvalid(PetscInt k);
  virtual void markColInvalid(PetscInt k);

  void localToGlobal(PetscInt k, PetscInt *i, PetscInt *j);

  virtual void addLocalResidualBlock(const PISMVector2 *y, PISMVector2 **yg);
  virtual void addLocalResidualBlock(const PetscScalar *y, PetscScalar **yg);

  virtual PetscErrorCode addLocalJacobianBlock(const PetscReal *K, Mat J);
  virtual PetscErrorCode setJacobianDiag(PetscInt i, PetscInt j, const PetscReal *K, Mat J);

  static const PetscInt Nk = 4; //<! The number of test functions defined on an element.
  
protected:
  static const PetscInt kDofInvalid = PETSC_MIN_INT / 8; //!< Constant for marking invalid row/columns.
  static const PetscInt kIOffset[Nk];
  static const PetscInt kJOffset[Nk];

  //! Indices of the current element (for updating residual/Jacobian).
  PetscInt m_i, m_j;

  //! Stencils for updating entries of the Jacobian matrix.
  MatStencil m_row[Nk], m_col[Nk]; 
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
with periodic boundaries can be thought of as using a kind of ghost.  So the rule still works:
compute on all elements containg a real vertex, but only update rows corresponding to that real vertex.

The calculation of what elements to index over needs to account for ghosts and the
presence or absense of periodic boundary conditions in the IceGrid.  The FEElementMap performs
that computation for you (see FEElementMap::xs and friends).

See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
*/
class FEElementMap
{
public:
  FEElementMap(const IceGrid &g);
  
  /*!\brief The total number of elements to be iterated over.  Useful for creating per-element storage.*/
  PetscInt element_count()
  {
    return xm*ym;
  }
  
  /*!\brief Convert an element index (\a i,\a j) into a flattened (1-d) array index, with the first
      element (\a i, \a j) to be iterated over corresponding to flattened index 0. */
  PetscInt flatten(PetscInt i, PetscInt j)
  {
    return (i-xs)*ym+(j-ys);
  }
  
  bool is_local_element(PetscInt i, PetscInt j, const IceGrid &g);
  
  PetscInt xs, //!< x-coordinate of the first element to loop over.
           xm, //!< total number of elements to loop over in the x-direction.
           ys, //!< y-coordinate of the first element to loop over.
           ym, //!< total number of elements to loop over in the y-direction. 
           lxs, //!< x-index of the first local element.
           lxm, //!< total number local elements in x direction.
           lys, //!< y-index of the first local element.
           lym; //!< total number local elements in y direction.

};

//! Numerical integration of finite element functions.
/*! The core of the finite element method is the computation of integrals over elements.
For nonlinear problems, or problems with non-constant coefficients (%i.e. any real problem)
the integration has to be done approximately:
\f[
\int_E f(x)\; dx \approx \sum_q f(x_q) w_q
\f]
for certain quadrature points \f$x_q\f$ and weights \f$w_q\f$.  An FEQuadrature is used
to evaluate finite element functions at quadrature points, and to compute weights \f$w_q\f$
for a given element.

In this concrete implementation, the reference element \f$R\f$ is the square 
\f$[-1,1]\times[-1,1]\f$.  On a given element, nodes (o) and quadrature points (*) 
are ordered as follows:
\verbatim
         3 o------------------o  2
           |  3             2 |
           |    *        *    |
           |                  |
           |                  |
           |    *        *    |
           |  0            1  |
        0  o------------------o  1
\endverbatim
So there are four quad points per element, which occur at \f$x,y=\pm 1/\sqrt{3}\f$.  This corresponds to the tensor product
of Gaussian integration on an interval that is exact for cubic functions on the interval.

Integration on a physical element can be thought of as being done by change of variables.  The quadrature weights need
to be modified, and the FEQuadrature takes care of this for you.  Because all elements in an IceGrid are congruent, the
quadrature weights are the same for each element, and are computed upon initialization with an IceGrid.

See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
*/
class FEQuadrature
{
public:

  static const PetscInt Nq = 4;  //!< Number of quadrature points.
  static const PetscInt Nk = 4;  //!< Number of test functions on the element.
  
  FEQuadrature();

  void init(const IceGrid &g,PetscScalar L=1.); // FIXME Allow a length scale to be specified.

  const FEFunctionGerm (*testFunctionValues())[Nq];  
  const FEFunctionGerm *testFunctionValues(PetscInt q);
  const FEFunctionGerm *testFunctionValues(PetscInt q,PetscInt k);
  

  void computeTrialFunctionValues( const PetscReal *x, PetscReal *vals);
  void computeTrialFunctionValues( const PetscReal *x, PetscReal *vals, PetscReal *dx, PetscReal *dy);
  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PetscReal const*const*xg, PetscReal *vals);
  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PetscReal const*const*xg, 
                                   PetscReal *vals, PetscReal *dx, PetscReal *dy);

  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals);
  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PetscReal (*Dv)[3]);  
  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PISMVector2 *dx, PISMVector2 *dy);  
  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PISMVector2 const*const*xg,  
                                   PISMVector2 *vals);
  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PISMVector2 const*const*xg,  
                                   PISMVector2 *vals, PetscReal (*Dv)[3]);

  void getWeightedJacobian(PetscReal *jxw);

  //! The coordinates of the quadrature points on the reference element.
  static const PetscReal quadPoints[Nq][2];
  //! The weights for quadrature on the reference element.
  static const PetscReal quadWeights[Nq];

protected:

  //! The Jacobian determinant of the map from the reference element to the physical element.
  PetscReal m_jacobianDet;
  //! Shape function values (for each of \a Nq quadrature points, and each of \a Nk shape function )
  FEFunctionGerm m_germs[Nq][Nk];

  PetscReal   m_tmpScalar[Nk];
  PISMVector2 m_tmpVector[Nk];
};


#endif/* _FETOOLS_H_*/
