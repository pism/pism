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

// Utility functions used by the SSAFEM code.

#include "SSAFEM_util.hh"


//  Integrations with integrands involving derivatives of finite element functions need a little extra work.  Let $f$ be a 
//  finite-element function (f=\sum_{ij} c_{ij} \psi_{ij}).  Consider a point p in physical element E, so p=F(q) where 
//  q is a point in the reference element and F is the map from reference to physical element.  Let g=f\circ F, so
//  
//       g = \sum_{n=1}^4 d_n \phi_n 
//
//  Then f(p) = f(F(q)) = g(q) = g(F^{-1}(p)) = \sum_{n=1}^4 d_n \phi_n(F^{-1}(p)).  The derivative of f with respect to $X$ or $Y$
//  can then be written in terms of derivatives of reference basis elements \phi_n and the Jacobian matrix of F^{-1} using the 
//  chain rule.  In particular, the derivative of f with respect to X at the point F(q_i) is
//
//     d f                            d\phi_n      d F^{-1}
//     --- (F(q_i))= \sum_{n=1}^4 d_n -----(q_i)  --------
//     d X                              dx           dX
//
//  ( we have used the fact that dF^{-1}/dY = 0, otherwise there would be an extra term). 
//
//  Similarly, 
//     d f                            d\phi_n      d F^{-1}
//     --- (F(q_i))= \sum_{n=1}^4 d_n -----(q_i)  --------
//     d Y                              dy           dY
//
//  In the code, the numbers dF^{-1}/dX and dF^{-1}/dY are referred to as the array jinvDiag.


FEElementMap::FEElementMap(const IceGrid &g)
{
  // Start by assuming ghost elements exist in all directions.
  // If there are 'mx' by 'my' points in the grid and no 
  // ghosts, there are (mx-1) by (my-1) elements.  In the presence
  // of ghosts there are then (mx+1) by (my+1) elements.
  xs= g.xs-1; xm = g.xm+1;
  ys= g.ys-1; ym = g.ym+1;

  // Now correct if needed. The only way there will not be ghosts is if the 
  // grid is not periodic and we are up against the grid boundary.
  
  if( !(g.periodicity & X_PERIODIC) )
  {
    // First element has x-index 0.
    if(xs < 0){
      xs = 0;
    }
    // Total number of elements is Mx-1, so xs+xm should be no larger.
    if(xs+xm > g.Mx-1) {
      xm = g.Mx-1-xs;
    }
  }

  if( !(g.periodicity & Y_PERIODIC) )
  {
    if(ys < 0){
      ys = 0;
    }
    if(ys+ym > g.My-1) {
      ym = g.My-1-ys;
    }
  }
}



//  Integrations with integrands involving derivatives of finite element functions need a little extra work.  Let $f$ be a 
//  finite-element function (f=\sum_{ij} c_{ij} \psi_{ij}).  Consider a point p in physical element E, so p=F(q) where 
//  q is a point in the reference element and F is the map from reference to physical element.  Let g=f\circ F, so
//  
//       g = \sum_{n=1}^4 d_n \phi_n 
//
//  Then f(p) = f(F(q)) = g(q) = g(F^{-1}(p)) = \sum_{n=1}^4 d_n \phi_n(F^{-1}(p)).  The derivative of f with respect to $X$ or $Y$
//  can then be written in terms of derivatives of reference basis elements \phi_n and the Jacobian matrix of F^{-1} using the 
//  chain rule.  In particular, the derivative of f with respect to X at the point F(q_i) is
//
//     d f                            d\phi_n      d F^{-1}
//     --- (F(q_i))= \sum_{n=1}^4 d_n -----(q_i)  --------
//     d X                              dx           dX
//
//  ( we have used the fact that dF^{-1}/dY = 0, otherwise there would be an extra term). 
//
//  Similarly, 
//     d f                            d\phi_n      d F^{-1}
//     --- (F(q_i))= \sum_{n=1}^4 d_n -----(q_i)  --------
//     d Y                              dy           dY
//
//  In the code, the numbers dF^{-1}/dX and dF^{-1}/dY are referred to as the array jinvDiag.



/*! Extract local degrees of freedom for element (\a i,\a j) from global vector \a xg to
  local vector \a x (scalar-valued DOF version). */
void FEDOFMap::extractLocalDOFs(PetscInt i,PetscInt j, PetscReal const*const*xg,PetscReal *x) const
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];  
}
/*! Extract local degrees of freedom for element (\a i,\a j) from global vector \a xg to
local vector \a x (vector-valued DOF version).
*/
void FEDOFMap::extractLocalDOFs(PetscInt i,PetscInt j, PISMVector2 const*const*xg,PISMVector2 *x) const
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];    
}
/*! Extract local degrees of freedom corresponding to the element set previously with \a reset. 
(scalar version)
*/
void FEDOFMap::extractLocalDOFs(PetscReal const*const*xg,PetscReal *x) const
{
  extractLocalDOFs(m_i,m_j,xg,x);
}
//! Extract local degrees of freedom for the element set previously with \a reset. (vector version)
void FEDOFMap::extractLocalDOFs(PISMVector2 const*const*xg,PISMVector2 *x) const
{
  extractLocalDOFs(m_i,m_j,xg,x);
}

//! Convert a local degree of freedom index \a k to a global degree of freedom index (\a i,\a j).
void FEDOFMap::localToGlobal(PetscInt k, PetscInt *i, PetscInt *j)
{
  *i = m_i + kIOffset[k];
  *j = m_j + kJOffset[k];  
}

/*! Initialize the FEDOFMap to element (\a i, \a j) for the purposes of inserting into
global residual and Jacobian arrays. */
void FEDOFMap::reset(PetscInt i, PetscInt j, const IceGrid &grid)
{
  m_i = i; m_j = j;
  // The meaning of i and j for a PISM IceGrid and for a Petsc DA are swapped (the so-called
  // fundamental transpose.  The interface between PISM and Petsc is the stencils, so all
  // interactions with the stencils involve a transpose.
  m_col[0].j = i;   m_col[0].i = j;
  m_col[1].j = i+1; m_col[1].i = j;
  m_col[2].j = i+1; m_col[2].i = j+1;
  m_col[3].j = i;   m_col[3].i = j+1;

  memcpy(m_row,m_col,Nk*sizeof(m_col[0]));

  // We do not ever sum into rows that are not owned by the local rank.
  for(PetscInt k=0; k<Nk; k++) {         
    PetscInt pism_i = m_row[k].j, pism_j = m_row[k].i;
    if (  pism_i < grid.xs || grid.xs+grid.xm-1 < pism_i || 
                pism_j < grid.ys || grid.ys+grid.ym-1 < pism_j ) {
      markRowInvalid(k);      
    }
  }
}

/*! Mark that the row corresponding to local degree of freedom \a k should not be updated
when inserting into the global residual or Jacobian arrays. */
void FEDOFMap::markRowInvalid(PetscInt k)
{
  m_row[k].i=m_row[k].j = kDofInvalid;
}

/*! Mark that the column corresponding to local degree of freedom \a k should not be updated
when inserting into the global Jacobian arrays. */
void FEDOFMap::markColInvalid(PetscInt k)
{
  m_col[k].i=m_col[k].j = kDofInvalid;
}

/*! Add the values of element-local residual contributions \a y to the global residual
 vector \a yg. */
void FEDOFMap::addLocalResidualBlock(const PISMVector2 *y, PISMVector2 **yg)
{
  for (int k=0; k<Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    yg[m_row[k].j][m_row[k].i].u += y[k].u;
    yg[m_row[k].j][m_row[k].i].v += y[k].v;
  }
}

//! Add the contributions of an element-local jaobian to the global Jacobian vector.
/*! The element-local Jacobian should be an array of Nk*Nk values in the
scalar case or (2Nk)*(2Nk) values in the vector valued case. */
PetscErrorCode FEDOFMap::addLocalJacobianBlock(const PetscReal *K, Mat J)
{
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,Nk,m_row,Nk,m_col,K,ADD_VALUES);CHKERRQ(ierr);  
  return 0;
}

//! Set a diagonal entry for global degree of freedom (\a i ,\a j ) in a Jacobian matrix
/* This is an unhappy hack for supporting Dirichlet constrained degrees of freedom.
In the scalar valued case, \a K should point to a single value, and in the vector case,
it should point to 4 (=2x2) values for the (2x2) block correspoinding to to u-u, u-v, v-u, and v-v
interactions at grid point (\a i, \aj).  Sheesh.*/
PetscErrorCode FEDOFMap::setJacobianDiag(PetscInt i, PetscInt j, const PetscReal*K, Mat J)
{
  MatStencil row;
  row.i=j; row.j=i;
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,1,&row,1,&row,K,INSERT_VALUES);CHKERRQ(ierr);  
  return 0;
}


FEQuadrature::FEQuadrature()
{  
}

void FEQuadrature::getWeightedJacobian(PetscReal *jxw)
{
  for(int q=0;q<Nq;q++)
  {
    jxw[q] = m_jacobianDet * quadWeights[q];
  } 
}

void FEQuadrature::init(const IceGrid &grid)
{
  // Since we use uniform cartesian coordinates, the Jacobian is constant and diagonal on every element.
  // Note that the reference element is \f$ [-1,1]^2 \f$ hence the extra factor of 1/2.
  PetscReal jacobian_x = 0.5*grid.dx;///ref.Length();
  PetscReal jacobian_y = 0.5*grid.dy;///ref.Length();
  m_jacobianDet = jacobian_x*jacobian_y;

  for(int q=0; q<Nq; q++){
    for(int k=0; k<Nk; k++){
      PetscInt qk = q*Nk+k;
      m_germs[q][k].val = interp[qk];
      m_germs[q][k].dx = derivx[qk]/jacobian_x;
      m_germs[q][k].dy = derivy[qk]/jacobian_y;
    }
  }
}

const FEFunctionGerm (*FEQuadrature::testFunctionValues())[4]
{
  return m_germs;
}

const FEFunctionGerm *FEQuadrature::testFunctionValues(PetscInt q)
{
  return m_germs[q];
}

const FEFunctionGerm *FEQuadrature::testFunctionValues(PetscInt q, PetscInt k)
{
  return m_germs[q] + k;
}

void FEQuadrature::computeTrialFunctionValues(const PetscReal *x, PetscReal *vals, PetscReal *dx, PetscReal *dy)
{
  for (int q=0; q<Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0; dx[q] = 0; dy[q] = 0;
    for (int k=0; k<Nk; k++) {
      vals[q] += test[k].val * x[k];
      dx[q]   += test[k].dx * x[k];
      dy[q]   += test[k].dy * x[k];
    }
  }
}
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PetscReal const*const*xg, 
                                 PetscReal *vals, PetscReal *dx, PetscReal *dy)
{
  dof.extractLocalDOFs(i,j,xg,m_tmpScalar);
  computeTrialFunctionValues(m_tmpScalar,vals,dx,dy);  
}


void FEQuadrature::computeTrialFunctionValues(const PetscReal *x, PetscReal *vals)
{
  for (int q=0; q<Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0;
    for (int k=0; k<Nk; k++) {
      vals[q] += test[k].val * x[k];
    }
  }
}
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, 
                                             PetscReal const*const*xg, PetscReal *vals)
{
  dof.extractLocalDOFs(i,j,xg,m_tmpScalar);
  computeTrialFunctionValues(m_tmpScalar,vals);
}


void FEQuadrature::computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PetscReal (*Dv)[3] )
{
  for (int q=0; q<Nq; q++) {
    vals[q].u = 0; vals[q].v = 0;
    PetscReal *Dvq = Dv[q];
    Dvq[0]=0; Dvq[1]=0; Dvq[2]=0;
    const FEFunctionGerm *test = m_germs[q];
    for(int k=0; k<Nk; k++) {
      vals[q].u += test[k].val * x[k].u;
      vals[q].v += test[k].val * x[k].v;
      Dvq[0] += test[k].dx * x[k].u;
      Dvq[1] += test[k].dy * x[k].v;
      Dvq[2] += 0.5*(test[k].dy*x[k].u + test[k].dx*x[k].v);
    }
  }  
}
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof,                                              
                                             PISMVector2 const*const*xg, PISMVector2 *vals, PetscReal (*Dv)[3] )
{
  dof.extractLocalDOFs(i,j,xg,m_tmpVector);
  computeTrialFunctionValues(m_tmpVector,vals,Dv);
}


void FEQuadrature::computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals)
{
  for (int q=0; q<Nq; q++) {
    vals[q].u = 0; vals[q].v = 0;
    const FEFunctionGerm *test = m_germs[q];
    for (int k=0; k<Nk; k++) {
      vals[q].u += test[k].val * x[k].u;
      vals[q].v += test[k].val * x[k].v;
    }
  }  
}
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof,                                              
                                             PISMVector2 const*const*xg, PISMVector2 *vals )
{
  dof.extractLocalDOFs(i,j,xg,m_tmpVector);
  computeTrialFunctionValues(m_tmpVector,vals);
}


//! Legacy code that needs to vanish. \todo Make it go away.
PetscTruth Floating(const IceFlowLaw &ice, PetscScalar ocean_rho,
                           PetscReal H, PetscReal bed)
{
  return ice.rho*H + ocean_rho*bed < 0 ? PETSC_TRUE : PETSC_FALSE;
}


//! Legacy code that needs to vanish. \todo Make it go away.
int PismIntMask(PetscScalar maskvalue) {
  return static_cast<int>(floor(maskvalue + 0.5));
}

