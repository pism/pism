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

PetscTruth Floating(const IceFlowLaw &ice, PetscScalar ocean_rho,
                           PetscReal H, PetscReal bed)
{
  return ice.rho*H + ocean_rho*bed < 0 ? PETSC_TRUE : PETSC_FALSE;
}


int PismIntMask(PetscScalar maskvalue) {
  return static_cast<int>(floor(maskvalue + 0.5));
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


IceGridElementIndexer::IceGridElementIndexer(const IceGrid &g)
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

// PetscErrorCode QuadGetStencils(DALocalInfo *info,PetscInt i,PetscInt j,
//                                       MatStencil row[],MatStencil col[])
// {
//   PetscErrorCode ierr;
//   PetscInt k;
// 
//   PetscFunctionBegin;
//   // PETSc reverses the meaning of i and j in stencil and info objects
//   ierr = PetscMemzero(col,4*sizeof(col[0]));CHKERRQ(ierr); // Probably just paranoia
//   col[0].j = i;   col[0].i = j;
//   col[1].j = i+1; col[1].i = j;
//   col[2].j = i+1; col[2].i = j+1;
//   col[3].j = i;   col[3].i = j+1;
//   ierr = PetscMemcpy(row,col,4*sizeof(col[0]));CHKERRQ(ierr);
//   for (k=0; k<4; k++) {         // We do not sum into rows that are not owned by us
//     if (   row[k].i < info->xs || info->xs+info->xm-1 < row[k].i
//         || row[k].j < info->ys || info->ys+info->ym-1 < row[k].j) {
//           // FIXME (DAM 2/11/11): Find the right negative number to use to indicate an invalid row/column.  
//           row[k].j = row[k].i = PETSC_MIN_INT/10;
//     }
//   }
//   PetscFunctionReturn(0);
// }


void DOFMap::extractLocalDOFs(PetscInt i,PetscInt j, PetscReal const*const*xg,PetscReal *x) const
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];  
}
void DOFMap::extractLocalDOFs(PetscInt i,PetscInt j, PISMVector2 const*const*xg,PISMVector2 *x) const
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];    
}

void DOFMap::extractLocalDOFs(PetscReal const*const*xg,PetscReal *x) const
{
  extractLocalDOFs(m_i,m_j,xg,x);
}
void DOFMap::extractLocalDOFs(PISMVector2 const*const*xg,PISMVector2 *x) const
{
  extractLocalDOFs(m_i,m_j,xg,x);
}

void DOFMap::reset(PetscInt i, PetscInt j, const IceGrid &grid)
{
  m_i = i; m_j = j;
  m_col[0].j = i;   m_col[0].i = j;
  m_col[1].j = i+1; m_col[1].i = j;
  m_col[2].j = i+1; m_col[2].i = j+1;
  m_col[3].j = i;   m_col[3].i = j+1;

  memcpy(m_row,m_col,Nk*sizeof(m_col[0]));

  for(PetscInt k=0; k<4; k++) {         // We do not sum into rows that are not owned by us
    PetscInt pism_i = m_row[k].j, pism_j = m_row[k].i;
    if (  pism_i < grid.xs || grid.xs+grid.xm-1 < pism_i || 
                pism_j < grid.ys || grid.ys+grid.ym-1 < pism_j ) {
      markRowInvalid(k);      
    }
  }
}

void DOFMap::localToGlobal(PetscInt k, PetscInt *i, PetscInt *j)
{
  *i = m_i + kIOffset[k];
  *j = m_j + kJOffset[k];  
}

void DOFMap::markRowInvalid(PetscInt k)
{
  m_row[k].i=m_row[k].j = kDofInvalid;
}
void DOFMap::markColInvalid(PetscInt k)
{
  m_col[k].i=m_col[k].j = kDofInvalid;
}

void DOFMap::addGlobalDOFs(const PISMVector2 *y, PISMVector2 **yg)
{
  for (int k=0; k<Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    yg[m_row[k].j][m_row[k].i].u += y[k].u;
    yg[m_row[k].j][m_row[k].i].v += y[k].v;
  }
}

PetscErrorCode DOFMap::addInteractionMatrix(const PetscReal *K, Mat J)
{
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,Nk,m_row,Nk,m_col,K,ADD_VALUES);CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode DOFMap::setDiag(PetscInt i, PetscInt j, const PetscReal*K, Mat J)
{
  MatStencil row;
  row.i=j; row.j=i;
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,1,&row,1,&row,K,INSERT_VALUES);CHKERRQ(ierr);  
  return 0;
}


Quadrature::Quadrature()
{  
}

void Quadrature::getWeightedJacobian(PetscReal *jxw)
{
  for(int q=0;q<Nq;q++)
  {
    jxw[q] = m_jacobianDet * quadWeights[q];
  } 
}

void Quadrature::init(const IceGrid &grid)
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

const FEFunctionGerm (*Quadrature::testFunctionValues())[4]
{
  return m_germs;
}

const FEFunctionGerm *Quadrature::testFunctionValues(PetscInt q)
{
  return m_germs[q];
}

const FEFunctionGerm *Quadrature::testFunctionValues(PetscInt q, PetscInt k)
{
  return m_germs[q] + k;
}

void Quadrature::computeTrialFunctionValues(const PetscReal *x, PetscReal *vals, PetscReal *dx, PetscReal *dy)
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
void Quadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof, PetscReal const*const*xg, 
                                 PetscReal *vals, PetscReal *dx, PetscReal *dy)
{
  dof.extractLocalDOFs(i,j,xg,m_tmpScalar);
  computeTrialFunctionValues(m_tmpScalar,vals,dx,dy);  
}


void Quadrature::computeTrialFunctionValues(const PetscReal *x, PetscReal *vals)
{
  for (int q=0; q<Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0;
    for (int k=0; k<Nk; k++) {
      vals[q] += test[k].val * x[k];
    }
  }
}
void Quadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof, 
                                             PetscReal const*const*xg, PetscReal *vals)
{
  dof.extractLocalDOFs(i,j,xg,m_tmpScalar);
  computeTrialFunctionValues(m_tmpScalar,vals);
}


void Quadrature::computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PetscReal (*Dv)[3] )
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
void Quadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof,                                              
                                             PISMVector2 const*const*xg, PISMVector2 *vals, PetscReal (*Dv)[3] )
{
  dof.extractLocalDOFs(i,j,xg,m_tmpVector);
  computeTrialFunctionValues(m_tmpVector,vals,Dv);
}


void Quadrature::computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals)
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
void Quadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const DOFMap &dof,                                              
                                             PISMVector2 const*const*xg, PISMVector2 *vals )
{
  dof.extractLocalDOFs(i,j,xg,m_tmpVector);
  computeTrialFunctionValues(m_tmpVector,vals);
}
