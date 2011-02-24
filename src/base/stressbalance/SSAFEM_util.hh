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

#ifndef _SSAFEM_UTIL_H_
#define _SSAFEM_UTIL_H_


#include <petscmat.h>
#include "iceModelVec.hh"
#include "flowlaws.hh"

//! \file 
//! \brief Classes for implementing the Finite Element Method on an IceGrid.
/*! \file
We assume that the reader has a basic understanding of the finite element method
at the level of ?????.  The following is a reminder of the method that also
gives the background for the how to implement it on an IceGrid with the tools in this 
module.

The IceGrid domain \Omega is decomposed into a grid of rectangular physical elements indexed by indices (i,j):

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

1) Evaluation of the residuals \f$r_{ij} = G_h(w_h,\psi_{ij})\f$

2) Evaluation of the Jacobian matrix
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

\li Extract from the global degrees of freedom \f$c\f$ defining \f$w_h\f$
the local degrees of freedom \f$d\f$ defining \f$\hat w_h\f$. (FEDOFMap)
\li Evaluate the local trial function \f$w_h\f$ (values and derivatives as needed)
at the quadrature points \f$x_q\f$ (FEQuadrature)
\li Evaluate the local test functions (again values and derivatives)
at the quadrature points. (FEQuadrature)
\li Obtain the quadrature weights $j_q$ for the element (FEQuadrature)
\liCompute the values of the integrand \f$g(\hat w_h,\psi_k)\f$
  at each quadrature point (call these \f$g_{qk}\f$) and
  form the weighted sums \f$y_k=\sum_{q} j_q g_{qk}\f$.
\li Each sum \f$y_k\f$ is the contribution of the current element to
a residual entry \f$r_{ij}\f$, where local degree of freedom \f$k\f$
corresponds with global degree of freedom \f$(i,j)\f$. The
local contibutions now need to be added to the global residual vector (FEDOFMap).

Computation of the Jacobian is similar, except that there are now
multiple integrals per element (one for each local degree of freedom of
\f$\hat w_h\f$).

All of the above is a simplified description of what happens in practice.
The complications below treated by the following classes, and discussed
in their documentation:

\li Ghost elements (as well as periodic boundary conditions): FEElementMap
\li Dirichlet data: FEDOFMap
\li Vector valued functions: (FEDOFMap, FEQuadrature)

It should be mentioned that the classes in this module are not intended
to be a fancy finite element package.
Their purpose is to clarify the steps that occur in the computation of
residuals and Jacobians in SSAFEM, and to isolate and consolodate 
the hard steps so that they are not scattered about the code.
*/

//*****************************************************************************************************
//
//  The reference element R is the square [-1,1]x[-1,1].  On a given element, nodes (o) and quadrature points (*) 
//  are ordered as follows:
//
//
//          3 o------------------o  2
//            |  3             2 |
//            |    *        *    |
//            |                  |
// R ->       |                  |
//            |    *        *    |
//            |  0            1  |
//         0  o------------------o  1
//
// The indices of the nodes are also used for the local index of the basis functions that are supported on the element.  


// Integration of a function over the reference element \a R is performed approximately using Gaussian integration:
//
// \int_R f \approx \sum_{i=1}^n w_i f(p_q)
//
// where the w_i are weights and the p_q are quadrature points.  In this implementation we have four quad points per
// element. The quadrature points occur at x,y=\pm \sqrt{3}.  This corresponds to the tensor product
// of Gaussian integration on an interval [-1,1] with quad points x=\pm\sqrt{3}, which is exact for cubic functions on the interval.

static const PetscReal quadPoints[4][2] = {{ -0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573,  0.57735026918962573 },
                                           { -0.57735026918962573,  0.57735026918962573 }};

// The weights w_i for gaussian quadrature on the reference element with these quadrature points
static const PetscReal quadWeights[4]  = {1,1,1,1};


// There are four reference basis functions $\phi_i$ per element; function $i$ is equal to 1 at node $i$ and equals zero at the other nodes.
// In order to perform integration involving functions of the form \sum_{k=0}^3 c_i \phi_i we need the values and derivatives
// of the reference basis functions at the quadrature points.  These are tabulated below.

#undef H
#undef L
#undef M
#undef P
// A linear function equal to 1 at x=-1 and equal to 0 at x=1 is 'H' at -1/\sqrt{3} and 'L' at 1/\sqrt{3}
#define H 0.78867513459481287
#define L 0.21132486540518708
// A linear function equal to 1 at x=-1 and equal to 0 at x=1 has derivative M;
// A linear function equal to 0 at x=-1 and equal to 1 at x=1 has derivative P.
#define M (-0.5)
#define P (0.5)

// interp[4*q + n] gives the value of basis function \a n at quadrature point \a q
static const PetscReal interp[4*4] = {H*H,H*L,L*L,L*H,  L*H,H*H,H*L,L*L,  L*L,H*L,H*H,L*H,  H*L,L*L,L*H,H*H};

// derivx[4*q + n] gives the derivative in the \a x direction of basis function \a n at quad point \a q
static const PetscReal derivx[4*4] = {M*H,P*H,P*L,M*L,  M*H,P*H,P*L,M*L,  M*L,P*L,P*H,M*H,  M*L,P*L,P*H,M*H};

// derivy[4*q + n] gives the derivative in the \a y direction of basis function \a n at quad point \a q
static const PetscReal derivy[4*4] = {H*M,L*M,L*P,H*P,  L*M,H*M,H*P,L*P,  L*M,H*M,H*P,L*P,  H*M,L*M,L*P,H*P};

#undef H
#undef L
#undef M
#undef P


// Computes the closest integer to maskvalue, with integers of the form n/2 rounded up.
int PismIntMask(PetscScalar maskvalue);


static const PetscInt kIOffset[4] = { 0, 1, 1, 0 };
static const PetscInt kJOffset[4] = { 0, 0, 1, 1 };

//! The mapping from global to local degrees of freedom.
/*! Computations of residual and Jacobian entries in the finite element method are
done by iterating of the elements and adding the various contributions from each element.  
To do this, degrees of freedom from the global vector of unknowns must be remapped to
element-local degrees of freedom for the purposes of local computation, 
(and the results added back again to the global residual and Jacobian arrays).

An FEDOFMap mediates the transfer between element-local and global degrees of freedom.

The FEDOFMap works equally well for vector and scalar valued finite element functions; the
degrees of freedom for a vector-value function are stored as PISMVector2's rather than
as PetscReal's.

The FEDOFMap is also (perhaps awkwardly) overloaded to also mediate transfering locally
computed contributions to residual and Jacobian matricies to their global
counterparts.

See also: \link SSAFEM_util.hh finite element tools background.\endlink.
*/
class FEDOFMap
{
public:
  FEDOFMap() {};
  
  virtual void extractLocalDOFs(PetscInt i, PetscInt j, PetscReal const*const*xg, PetscReal *x) const;
  virtual void extractLocalDOFs(PetscInt i, PetscInt j, PISMVector2 const*const*xg, PISMVector2 *x) const;

  virtual void extractLocalDOFs(PetscReal const*const*xg, PetscReal *x) const;
  virtual void extractLocalDOFs(PISMVector2 const*const*xg, PISMVector2 *x) const;

  virtual void reset(PetscInt i, PetscInt j, const IceGrid &g);
  
  virtual void markRowInvalid(PetscInt k);
  virtual void markColInvalid(PetscInt k);

  void localToGlobal(PetscInt k, PetscInt *i, PetscInt *j);

  virtual void addLocalResidualBlock(const PISMVector2 *y, PISMVector2 **yg);
  virtual PetscErrorCode addLocalJacobianBlock(const PetscReal *K, Mat J);
  virtual PetscErrorCode setJacobianDiag(PetscInt i, PetscInt j, const PetscReal *K, Mat J);

  static const PetscInt Nk = 4; //<! The number of test functions defined on an element.
  
protected:
  static const PetscInt kDofInvalid = PETSC_MIN_INT / 8; //<! Constant for marking invalid row/columns.

  PetscInt m_i, m_j; //<! index of the current element.
  MatStencil m_row[Nk], m_col[Nk]; //<! row and column indices for the global degrees of freedom for this element.
};


class FEElementMap
{
public:
  FEElementMap(const IceGrid &g);
  
  PetscInt element_count()
  {
    return (xm-xs)*(ym-ys);
  }
  
  PetscInt flatten(PetscInt i, PetscInt j)
  {
    return (i-xs)*ym+(j-ys);
  }
  
  PetscInt xs, //<! Start o
           xm, 
           ys, 
           ym;
  
};

struct FEFunctionGerm
{
  PetscReal val, dx, dy;
};


class FEQuadrature
{
public:

  static const PetscInt Nq = 4;  // Number of quadrature points.
  static const PetscInt Nk = 4;  // Number of test functions on the element.
  
  FEQuadrature();

  void init(const IceGrid &g); // FIXME Allow a length scale to be specified.

  const FEFunctionGerm (*testFunctionValues())[4];  
  const FEFunctionGerm *testFunctionValues(PetscInt q);
  const FEFunctionGerm *testFunctionValues(PetscInt q,PetscInt k);
  
  
  void computeTrialFunctionValues( const PetscReal *x, PetscReal *vals);
  void computeTrialFunctionValues( const PetscReal *x, PetscReal *vals, PetscReal *dx, PetscReal *dy);

  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PetscReal const*const*xg, PetscReal *vals);
  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PetscReal const*const*xg, 
                                   PetscReal *vals, PetscReal *dx, PetscReal *dy);

  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PetscReal (*Dv)[3]);
  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals);

  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PISMVector2 const*const*xg,  
                                   PISMVector2 *vals, PetscReal (*Dv)[3]);

  void computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PISMVector2 const*const*xg,  
                                   PISMVector2 *vals);

  void getWeightedJacobian(PetscReal *jxw);

protected:
  PetscReal m_jacobianDet;
  FEFunctionGerm m_germs[Nq][Nk];
  PetscReal   m_tmpScalar[Nk];
  PISMVector2 m_tmpVector[Nk];
};



// Returns true if bed<0 and the mass of a column of ice of height H  is less than the mass of
// a column of (sea-)water of height |bed|.
PetscTruth Floating(const IceFlowLaw &ice, PetscScalar ocean_rho,
                          PetscReal H, PetscReal bed);

#endif/* _SSAFEM_UTIL_H_*/
