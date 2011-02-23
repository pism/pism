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

void QuadZeroScalar(PetscReal x[])
{
  PetscMemzero(x,4*sizeof(x[0]));
}

void QuadZeroVel(PISMVector2 x[])
{
  PetscMemzero(x,4*sizeof(x[0]));
}

void QuadExtractScalar(PetscInt i,PetscInt j,PetscReal **xg,PetscReal x[])
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];
}
void QuadInsertScalar(PetscInt i,PetscInt j,PetscReal x[],PetscReal **xg)
{
  xg[i][j] = x[0]; xg[i+1][j] = x[1]; xg[i+1][j+1] = x[2]; xg[i][j+1] = x[3];
}

void QuadExtractVel(PetscInt i,PetscInt j,const PISMVector2 **xg,PISMVector2 x[])
{
  x[0].u = xg[i][j].u;        x[0].v = xg[i][j].v;
  x[1].u = xg[i+1][j].u;      x[1].v = xg[i+1][j].v;
  x[2].u = xg[i+1][j+1].u;    x[2].v = xg[i+1][j+1].v;
  x[3].u = xg[i][j+1].u;      x[3].v = xg[i][j+1].v;
}

void QuadInsertVel(const MatStencil row[],const PISMVector2 x[],PISMVector2 **xg)
{
#if 0
  printf("El(%d,%d): inserting [(%g,%g) (%g,%g) (%g,%g) (%g,%g)]\n",
         i,j,x[0].u,x[0].v,x[1].u,x[1].v,x[2].u,x[2].v,x[3].u,x[3].v);
#endif
  for (int k=0; k<4; k++) {
    if (row[k].i < 0 || row[k].j < 0) continue;
    // The stencil is in PETSc ordering [j][i]
    xg[row[k].j][row[k].i].u += x[k].u;
    xg[row[k].j][row[k].i].v += x[k].v;
  }
}

void QuadMatMultScalar(const PetscReal *A,const PetscReal *x,PetscReal *y)
{
  PetscInt i,j;
  const PetscReal *a;

  for (i=0,a=A; i<4; i++,a+=4) {
    y[i] = 0;
    for (j=0; j<4; j++) {
      y[i] += a[j] * x[j];
    }
  }
}

void QuadMatMultVel(const PetscReal *A,const PISMVector2 *x,PISMVector2 *y)
{
  PetscInt i,j;
  const PetscReal *a;

  for (i=0,a=A; i<4; i++,a+=4) {
    y[i].u = y[i].v = 0;
    for (j=0; j<4; j++) {
      y[i].u += a[j] * x[j].u;
      y[i].v += a[j] * x[j].v;
    }
  }
}

void QuadMatMultTransposeScalar(const PetscReal *A,const PetscReal *x,PetscReal *y)
{
  PetscInt i,j;
  const PetscReal *a;

  QuadZeroScalar(y);
  for (i=0,a=A; i<4; i++,a+=4) {
    for (j=0; j<4; j++) {
      y[j] += a[j] * x[i];
    }
  }
}

void QuadMatMultTransposeVel(const PetscReal *A,const PISMVector2 *x,PISMVector2 *y)
{
  PetscInt i,j;
  const PetscReal *a;

  QuadZeroVel(y);
  for (i=0,a=A; i<4; i++,a+=4) {
    for (j=0; j<4; j++) {
      y[j].u += a[j] * x[i].u;
      y[j].v += a[j] * x[i].v;
    }
  }
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

#undef __FUNCT__
#define __FUNCT__ "QuadEvaluateVel"
PetscErrorCode QuadEvaluateVel(const PISMVector2 *x,PetscInt q,
                                      const PetscReal jinvDiag[],
                                      PISMVector2 *u,PetscReal Du[])
{

  PetscFunctionBegin;
  u->u = u->v = 0;
  Du[0] = Du[1] = Du[2] = 0;
  for (PetscInt k=0; k<4; k++) {   // loop over trial functions to compute function values at quadrature point
    const PetscInt qk = q*4+k;
    u->u  += interp[qk] * x[k].u; u->v  += interp[qk] * x[k].v; // velocity at q
    Du[0] += jinvDiag[0]*derivx[qk] * x[k].u;                 // u_x
    Du[1] += jinvDiag[1]*derivy[qk] * x[k].v;                 // v_y
    Du[2] += 0.5*(jinvDiag[1]*derivy[qk] * x[k].u + jinvDiag[0]*derivx[qk] * x[k].v); // (u_y + v_x)/2
  }
  // FIXME: Turn these back on?
  // PismValidVelocity(*u);
  // PismValidStrainRate(Du);
  PetscFunctionReturn(0);
}

PetscErrorCode QuadGetStencils(DALocalInfo *info,PetscInt i,PetscInt j,
                                      MatStencil row[],MatStencil col[])
{
  PetscErrorCode ierr;
  PetscInt k;

  PetscFunctionBegin;
  // PETSc reverses the meaning of i and j in stencil and info objects
  ierr = PetscMemzero(col,4*sizeof(col[0]));CHKERRQ(ierr); // Probably just paranoia
  col[0].j = i;   col[0].i = j;
  col[1].j = i+1; col[1].i = j;
  col[2].j = i+1; col[2].i = j+1;
  col[3].j = i;   col[3].i = j+1;
  ierr = PetscMemcpy(row,col,4*sizeof(col[0]));CHKERRQ(ierr);
  for (k=0; k<4; k++) {         // We do not sum into rows that are not owned by us
    if (   row[k].i < info->xs || info->xs+info->xm-1 < row[k].i
        || row[k].j < info->ys || info->ys+info->ym-1 < row[k].j) {
          // FIXME (DAM 2/11/11): Find the right negative number to use to indicate an invalid row/column.  
          row[k].j = row[k].i = PETSC_MIN_INT/10;
    }
  }
  PetscFunctionReturn(0);
}
