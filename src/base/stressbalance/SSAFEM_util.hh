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



//*****************************************************************************************************
//
//  This module contains helper routines for implementing the finite element method
//  on a structured grid.  
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

static const PetscInt numQuadPoints = 4;
static const PetscReal quadPoints[4][2] = {{ -0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573,  0.57735026918962573 },
                                           { -0.57735026918962573,  0.57735026918962573 }};

// The weights w_i for gaussian quadrature on the reference element with these quadrature points
static const PetscReal quadWeights[4]  = {1,1,1,1};


// There are four reference basis functions $\phi_i$ per element; function $i$ is equal to 1 at node $i$ and equals zero at the other nodes.
// In order to perform integration involving functions of the form \sum_{k=0}^3 c_i \phi_i we need the values and derivatives
// of the reference basis functions at the quadrature points.  These are tabulated below.

// UNUSED:
//static const PetscReal quadWeights1[2] = {1,1};

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


// The global rectangular domain \Omega is decomposed into a grid of rectangular physical elements indexed by indices (i,j):
//
//      (0,1) (1,1)
//     ---------
//     |   |   |
//     ---------
//     |   |   |
//     ---------
//    (0,0) (1,0)   
//
//
// For each vertex (i,j) there is a basis function \psi_{ij} that is bilinear on each physical element, equals 1 at
// vertex (i,j), and equals 0 at all other vertices.  A finite element function \a f on the domain is written as a sum
//
// f = \sum_{i,j} c_{ij} \psi_{ij}

// For each physical element E_{ij}, there is an affine map F_{ij} from the reference element R to E_{ij}.  In this implementation,
// the rectangles in the domain are all congruent and the maps F_{ij} are all the same up to a translation:
//
//      F_{ij}(x,y) = ( dX*(x+1)/2, dY*(y+1)/2 ) + V_{ij} where dX and dY are the spacing between elements in the physical domain

//  and V_{ij} is the lower-left vertex of the physical element.  If \a f is a finite-element function on the whole domain,
//  then f_{ij} = f\circ F_{ij} is a function defined on the reference element and is a sum of element-basis functions.  

// There is a map from the degrees of freedom of a finite element function $f$ defined on the whole domain to the degrees of
// freedom of f_{ij}, its 'restriction' to the reference element. Specifically, element E_{ij} has vertices with global indices
//  (i,j), (i+1,j), (i+1,j+1), and (i,j+1).  These are mapped to reference degreees of freedom 0-3 in this order.

//  Integration on the whole domain is acheived by summing integrals over each physical element.  Integration on a physcial
//  element is converted to integration over the reference element, and integration on the reference element is done approximately
//  using Gaussian integration as described previously.  The conversion of an integral on a physical element to the reference
//  element proceeds as follows.  Let E be the physical element in question and F the map from R to E.  Let J be the Jacobian
//  matrix of F and let detJ be its determinant (i.e. dX*dY/4 for the specific case at hand).
//
//  Then    \int_E f(X,Y) dXdY   =   \int_R  f(F(x,y)) detJ dx dy  \approx \sum_{i} w_i detJ f(F(q_i))
//
//  So the key point in computing the integral is to obtain the values of f  at F(q_i), for each quadrature point.
//
//  For example, suppose that f is a finite element function on the entire domain (\Omega)  with global 
//  degrees of freedom c_{kl} and we wish to compute \int_\Omega f. This is a sum: 
//
//     \int_\Omega f(X,Y) dXdY = \sum_{ij} \int_{E_{ij}} f(X,Y) dXdY = \sum_{ij} \int_R  f(F_{ij}(x,y)) detJ dx dy.  
//  
//  To compute the approximate value of each term in the sum of the form \int_R f(F(x,y)) detJ dx dy:
//
//  1) Convert the global degree of freedom weights to element local degree of freedom weights: c_{kl} ->  d_n (n=0,1,2,3).  The
//     function f\circ F is then \sum_{n=0}^3 d_n \psi_n.
//
//  2) Compute the values of \sum_{n=0}^3 d_n \psi_n at each quadrature point q_i.  That is, compute the values
//     r_i = \sum_{n=0}^3 \psi_n(q_i).  Then f(F(q_i)) = r_i.
//
//  3) Compute the weighted sum over the integration domain: \sum_{i=1}^4 w_i detJ r_i
//
//  Step 1) is acheived with the function QuadExtractScalar.  
//  Step 2) is acheived using the function QuadMatMultScalar.  In particular, let A_{i,n} be the value of basis function \a n at 
//  quadrature point \a i.  Then r_i= \sum_n A_{i,n} d_n.  The matrix A_{i,n} is the variable \a interp above, and the matrix
//  multiplication is performed with QuadMatMultScalar.
//  Step 3) is simply done explicitly.
//
//  Integrations with integrands involving a more complication function of a finite-element function, or involving a combination
//  of several finite elements functions, are perfomed similarly.  Steps 1) and 2) need to be done for each of the finite-element
//  functions involved in the integrand, and the form of the sum in Step 3) needs to reflect the new integrand (the terms w_i and detJ
//  don't change, but the remainder will).  We can also compute integrands that involve derivatives of finite element functions;
//  the computation of the derivatives of a finite element function \a f at some F(q_i) is described in SSAFEM_util.cc.  The computation
//  is a little tedious, and is best not done explicitly.  See QuadEvaluateVel.
//



// // Given an element-local array of basis function weights \a x for element (i,j), adds the weights
// // into the global array of basis function weights \a xg. (x and xg are vectors)
// void QuadInsertVel(const MatStencil row[],const PISMVector2 x[],PISMVector2 **xg);
// 
// 
// 
// // Computes the closest integer to maskvalue, with integers of the form n/2 rounded up.
// int PismIntMask(PetscScalar maskvalue);
// 
// 
// PetscErrorCode QuadGetStencils(DALocalInfo *info,PetscInt i,PetscInt j,
//                                MatStencil row[],MatStencil col[]);


static const PetscInt kIOffset[4] = { 0, 1, 1, 0 };
static const PetscInt kJOffset[4] = { 0, 0, 1, 1 };

class DOFMap
{
public:
  DOFMap() {};
  virtual void extractLocalDOFs(PetscInt i, PetscInt j, PetscReal const*const*xg, PetscReal *x) const;
  virtual void extractLocalDOFs(PetscInt i, PetscInt j, PISMVector2 const*const*xg, PISMVector2 *x) const;

  virtual void extractLocalDOFs(PetscReal const*const*xg, PetscReal *x) const;
  virtual void extractLocalDOFs(PISMVector2 const*const*xg, PISMVector2 *x) const;

  virtual void reset(PetscInt i, PetscInt j, const IceGrid &g);
  
  virtual void markRowInvalid(PetscInt k);
  virtual void markColInvalid(PetscInt k);

  void localToGlobal(PetscInt k, PetscInt *i, PetscInt *j);

  virtual void addGlobalDOFs(const PISMVector2 *y, PISMVector2 **yg);
  virtual PetscErrorCode addInteractionMatrix(const PetscReal *K, Mat J);
  virtual PetscErrorCode setDiag(PetscInt i, PetscInt j, const PetscReal *K, Mat J);
  
protected:
  static const PetscInt Nk = 4;
  static const PetscInt kDofInvalid = PETSC_MIN_INT / 8;

  PetscInt m_i, m_j;
  MatStencil m_row[Nk];
  MatStencil m_col[Nk];
};


class IceGridElementIndexer 
{
public:
  IceGridElementIndexer(const IceGrid &g);
  
  PetscInt element_count()
  {
    return (xm-xs)*(ym-ys);
  }
  
  PetscInt flatten(PetscInt i, PetscInt j)
  {
    return (i-xs)*ym+(j-ys);
  }
  
  PetscInt xs, xm, ys, ym;
  
};

struct FEFunctionGerm
{
  PetscReal val, dx, dy;
};


class Quadrature
{
public:

  static const PetscInt Nq = 4;  // Number of quadrature points.
  static const PetscInt Nk = 4;  // Number of test functions on the element.
  
  Quadrature();

  void init(const IceGrid &g); // FIXME Allow a length scale to be specified.

  const FEFunctionGerm (*testFunctionValues())[4];  
  const FEFunctionGerm *testFunctionValues(PetscInt q);
  const FEFunctionGerm *testFunctionValues(PetscInt q,PetscInt k);
  
  
  void computeTrialFunctionValues( const PetscReal *x, PetscReal *vals);
  void computeTrialFunctionValues( const PetscReal *x, PetscReal *vals, PetscReal *dx, PetscReal *dy);

  void computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof, PetscReal const*const*xg, PetscReal *vals);
  void computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof, PetscReal const*const*xg, 
                                   PetscReal *vals, PetscReal *dx, PetscReal *dy);

  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PetscReal (*Dv)[3]);
  void computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals);

  void computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof, PISMVector2 const*const*xg,  
                                   PISMVector2 *vals, PetscReal (*Dv)[3]);

  void computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof, PISMVector2 const*const*xg,  
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
