// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cmath>
#include <cstring>
#include <petscda.h>
#include "iceModel.hh"


PetscErrorCode IceModel::setupForSSA(const PetscScalar minH) {
  PetscErrorCode ierr;
  PetscScalar **h, **H, **mask, **ubar, **vbar;
  PetscScalar ***u, ***v;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) != MASK_SHEET && H[i][j] < minH) {
        h[i][j] += (minH - H[i][j]);
        H[i][j] = minH;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceModel::cleanupAfterSSA(const PetscScalar minH) {
  PetscErrorCode  ierr;
  PetscScalar **h, **H, **mask;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) != MASK_SHEET && H[i][j] <= minH) {
        h[i][j] -= minH;
        H[i][j] = 0.0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  return 0;
}


//! Compute the product of the effective viscosity \f$\nu\f$ and ice thickness \f$H\f$.
/*! 
@cond CONTINUUM
In PISM the product \f$\nu H\f$ can be (i) constant, (ii) can be computed with a constant ice hardness
\f$\bar B\f$ (temperature-independent) but with dependence of the viscosity on the strain rates, or 
(iii) it can depend on the strain rates and have a vertically-averaged ice hardness.

The flow law in ice stream and ice shelf regions must, for now, be a temperature-dependent Glen
law.  (That is, more general forms like Goldsby-Kohlstedt are not yet inverted.)

The effective viscosity is
   \f[ \nu = \frac{\bar B}{2} \left[\left(\frac{\partial u}{\partial x}\right)^2 + 
                               \left(\frac{\partial v}{\partial y}\right)^2 + 
                               \frac{\partial u}{\partial x} \frac{\partial v}{\partial y} + 
                               \frac{1}{4} \left(\frac{\partial u}{\partial y}
                                                 + \frac{\partial v}{\partial x}\right)^2
                               \right]^{(1-n)/(2n)}                                                \f]
where in the temperature-dependent case
   \f[ \bar B = (\text{FIXME: SOME INTEGRAL})\f]
@endcond
@cond NUMERIC
In the temperature dependent case the ice hardness is computed by the trapezoid rule.
@endcond
 */
PetscErrorCode IceModel::computeEffectiveViscosity(Vec vNu[2], PetscReal epsilon) {
  PetscErrorCode ierr;

  if (leaveNuAloneSSA == PETSC_TRUE) {
    return 0;
  }
  if (useConstantNuForSSA == PETSC_TRUE) {
    ierr = VecSet(vNu[0], constantNuForSSA); CHKERRQ(ierr);
    ierr = VecSet(vNu[1], constantNuForSSA); CHKERRQ(ierr);
    return 0;
  }
  
  // this constant is the form of regularization used by C. Schoof 2006 "A variational
  // approach to ice streams" J Fluid Mech 556 pp 227--251
  const PetscReal  schoofReg = PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof);

  // We need to compute integrated effective viscosity. It is locally determined by the
  // strain rates and temperature field.
  PetscScalar **mask, **H, ***T, **nu[2], **u, **v;  
  ierr = VecSet(vNu[0], 0.0); CHKERRQ(ierr);
  ierr = VecSet(vNu[1], 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u);
  ierr = DAVecGetArray(grid.da2, vvbar, &v);
  for (PetscInt o=0; o<2; ++o) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscInt      oi = 1-o, oj=o;
        if (intMask(mask[i][j]) != MASK_SHEET || intMask(mask[i+oi][j+oj]) != MASK_SHEET) {
          const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
          PetscScalar u_x, u_y, v_x, v_y;
          // Check the offset to determine how to differentiate velocity
          if (o == 0) {
            u_x = (u[i+1][j] - u[i][j]) / dx;
            u_y = (u[i][j+1] + u[i+1][j+1] - u[i][j-1] - u[i+1][j-1]) / (4*dy);
            v_x = (v[i+1][j] - v[i][j]) / dx;
            v_y = (v[i][j+1] + v[i+1][j+1] - v[i][j-1] - v[i+1][j-1]) / (4*dy);
          } else {
            u_x = (u[i+1][j] + u[i+1][j+1] - u[i-1][j] - u[i-1][j+1]) / (4*dx);
            u_y = (u[i][j+1] - u[i][j]) / dy;
            v_x = (v[i+1][j] + v[i+1][j+1] - v[i-1][j] - v[i-1][j+1]) / (4*dx);
            v_y = (v[i][j+1] - v[i][j]) / dy;
          }

          const PetscScalar myH = 0.5 * (H[i][j] + H[i+oi][j+oj]);
          if (useConstantHardnessForSSA == PETSC_TRUE) { 
            // constant \bar B case, i.e for EISMINT ROSS; "nu" is really "nu H"!
            // FIXME: assumes n=3; can be fixed by having new ice.effectiveViscosityConstant() method
            nu[o][i][j] = myH * constantHardnessForSSA * 0.5 *
              pow(schoofReg + PetscSqr(u_x) + PetscSqr(v_y) + 0.25*PetscSqr(u_y+v_x) + u_x*v_y,
                  -(1.0/3.0));
          } else { 
            // usual temperature-dependent case; "nu" is really "nu H"!
            nu[o][i][j] = ice.effectiveViscosityColumn(schoofReg,
                                    myH,dz, u_x, u_y, v_x, v_y, T[i][j], T[i+oi][j+oj]);
          }

          if (! finite(nu[o][i][j]) || false) {
            ierr = PetscPrintf(grid.com, "nu[%d][%d][%d] = %e\n", o, i, j, nu[o][i][j]);
            CHKERRQ(ierr); 
            ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n", 
                               u_x, u_y, v_x, v_y);
            CHKERRQ(ierr);
          }
          
          // We ensure that nu is bounded below by a positive constant.
          nu[o][i][j] += epsilon;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr); CHKERRQ(ierr);

  // Some communication
  ierr = DALocalToLocalBegin(grid.da2, vNu[0], INSERT_VALUES, vNu[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vNu[0], INSERT_VALUES, vNu[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vNu[1], INSERT_VALUES, vNu[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vNu[1], INSERT_VALUES, vNu[1]); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::testConvergenceOfNu(Vec vNu[2], Vec vNuOld[2],
                                             PetscReal *norm, PetscReal *normChange) {
  PetscErrorCode  ierr;
  PetscReal nuNorm[2], nuChange[2];
  const PetscScalar area = grid.p->dx * grid.p->dy;
#define MY_NORM     NORM_2

  // Test for change in nu
  ierr = VecAXPY(vNuOld[0], -1, vNu[0]); CHKERRQ(ierr);
  ierr = VecAXPY(vNuOld[1], -1, vNu[1]); CHKERRQ(ierr);
  ierr = DALocalToGlobal(grid.da2, vNuOld[0], INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecNorm(g2, MY_NORM, &nuChange[0]); CHKERRQ(ierr);
  nuChange[0] *= area;
  ierr = DALocalToGlobal(grid.da2, vNuOld[1], INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecNorm(g2, MY_NORM, &nuChange[1]); CHKERRQ(ierr);
  nuChange[1] *= area;
  *normChange = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));

  ierr = DALocalToGlobal(grid.da2, vNu[0], INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecNorm(g2, MY_NORM, &nuNorm[0]); CHKERRQ(ierr);
  nuNorm[0] *= area;
  ierr = DALocalToGlobal(grid.da2, vNu[1], INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecNorm(g2, MY_NORM, &nuNorm[1]); CHKERRQ(ierr);
  nuNorm[1] *= area;
  *norm = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));
  return 0;
}


//! Assemble the matrix (left-hand side) for the numerical approximation of the SSA equations.
/*! 
@cond CONTINUUM
The SSA equations are in their clearest form
    \f[ - \frac{\partial T_{ij}}{\partial x_j} + \tau_{(b)i} = f_i \f]
where \f$i,j\f$ range over \f$x,y\f$, \f$T_{ij}\f$ is a depth-integrated viscous stress tensor 
(i.e. equation (2.6) in (Schoof 2006), and following (Morland 1987)), 
and \f$\tau_{(b)i}\f$ are the components of the basal shear stress.  
Also \f$f_i\f$ is the driving shear stress \f$f_i = - \rho g H \frac{\partial h}{\partial x_i}\f$.  
These equations determine velocity in a more-or-less elliptic equation manner.  Here \f$H\f$ 
is the ice thickness and \f$h\f$ is the elevation of the surface of the ice.

More specifically, the SSA equations are
    \f[ - 2 \frac{\partial}{\partial x}\left[\nu H \left(2 \frac{\partial u}{\partial x}
                                                       + \frac{\partial v}{\partial y}\right)\right]
        - \frac{\partial}{\partial y}\left[\nu H \left(\frac{\partial u}{\partial y}
                                                       + \frac{\partial v}{\partial x}\right)\right]
        + \tau_{(b)x}
        = - \rho g H \frac{\partial h}{\partial x}, \f]
    \f[ - \frac{\partial}{\partial x}\left[\nu H \left(\frac{\partial u}{\partial y}
                                                       + \frac{\partial v}{\partial x}\right)\right]
        - 2 \frac{\partial}{\partial y}\left[\nu H \left(\frac{\partial u}{\partial x}
                                                         + 2 \frac{\partial v}{\partial y}\right)\right]
        + \tau_{(b)y}
        = - \rho g H \frac{\partial h}{\partial y}, \f]
where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the \f$y\f$-component 
of the velocity.  Note \f$\nu\f$ is the vertically-averaged effective viscosity of the ice.  

For ice shelves \f$\tau_{(b)i} = 0\f$ (MacAyeal et al 1996).  For ice streams with a basal 
till modelled as a plastic material, \f$\tau_{(b)i} = \tau_c u_i/|\mathbf{u}|\f$ where 
\f$\mathbf{u} = (u,v)\f$, \f$|\mathbf{u}| = \left(u^2 + v^2\right)^{1/2}\f$, and \f$\tau_c\f$ 
is the yield stress of the till  (Schoof 2006).  For ice streams with a basal till modelled 
as a linearly-viscous material, \f$\tau_{(b)i} = \beta u_i\f$ where \f$\beta\f$ is the basal
drag (friction) parameter (Hulbe & MacAyeal 1999).

@endcond
@cond NUMERIC
Note that the basal shear stress appears on the \emph{left} side of the above system.  
We believe this is crucial, because of its effect on the spectrum of the linear 
approximations of each stage.  The effect on spectrum is clearest in the linearly-viscous
till case (i.e. Hulbe & MacAyeal 1999) but there seems to be an analogous effect in the 
plastic till case (Schoof 2006).

This method assembles the matrix for the left side of the SSA equations.  The numerical method is finite difference.  In particular [FIXME: explain f.d. approxs, esp. mixed derivatives]
@endcond
 */
PetscErrorCode IceModel::assembleSSAMatrix(Vec vNu[2], Mat A) {
  const PetscInt  Mx=grid.p->Mx, My=grid.p->My, M=2*My;
  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy;
  const PetscScalar   one = 1.0;
  PetscErrorCode  ierr;
  PetscScalar     **mask, **nu[2], **u, **v, **beta, **tauc;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  /* matrix assembly loop */
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j;
      const PetscInt rowU = i*M + J;
      const PetscInt rowV = i*M + J+1;
      if (intMask(mask[i][j]) == MASK_SHEET) {
        ierr = MatSetValues(A, 1, &rowU, 1, &rowU, &one, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, 1, &rowV, &one, INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscInt im = (i + Mx - 1) % Mx, ip = (i + 1) % Mx;
        const PetscInt Jm = 2 * ((j + My - 1) % My), Jp = 2 * ((j + 1) % My);
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up. */

// these DO use thickness H because nu[][][] is actually viscosity times thickness
        const PetscScalar c00 = nu[0][i-1][j];
        const PetscScalar c01 = nu[0][i][j];
        const PetscScalar c10 = nu[1][i][j-1];
        const PetscScalar c11 = nu[1][i][j];

#if 0
        /*
        * These are stencils for
        *  - \nabla \cdot ( \nu \nabla u ) = 0
        *      and
        *  - \nabla \cdot ( \nu \nabla v ) = 0
        * They leave u and v uncoupled.
        */
        const PetscInt stencilSize = 5;
        const PetscInt colU[stencilSize] = {
          i*M+Jp,
          im*M+J,     i*M+J,      ip*M+J,
          i*M+Jm};
        const PetscScalar valU[stencilSize] = {
          -c11/dy2,
          -c00/dx2,   (c00+c01)/dx2+(c10+c11)/dy2,    -c01/dx2,
          -c10/dy2};
        const PetscInt colV[stencilSize] = {
          i*M+Jp+1,
          im*M+J+1,   i*M+J+1,    ip*M+J+1,
          i*M+Jm+1};
        const PetscScalar valV[stencilSize] = {
          -c11/dy2,
          -c00/dx2,   (c00+c01)/dx2+(c10+c11)/dy2,    -c01/dx2,
          -c10/dy2};
#else
#if 0
        /*
        * These are the stencils for constant thickness and constant viscosity.
        * They are not scaled for grid size so use 0 as RHS.
        * Only use for basic testing.
        */
        const PetscInt stencilSize = 9;
        const PetscInt colU[] = {
          im*M+Jp+1,  i*M+Jp,     ip*M+Jp+1,
          im*M+J,     i*M+J,      ip*M+J,
          im*M+Jm+1,  i*M+Jm,     ip*M+Jm+1 };
        const PetscScalar valU[] = {0.75, -1, -0.75,        -4, 10, -4,     -0.75, -1, 0.75};
        const PetscInt colV[] = {
          im*M+Jp,    i*M+Jp+1,   ip*M+Jp,
          im*M+J+1,   i*M+J+1,    ip*M+J+1,
          im*M+Jm,    i*M+Jm+1,   ip*M+Jm };
        const PetscScalar valV[] = {0.75, -4, -0.75,        -1, 10, -1,     -0.75, -4, 0.75};
#else
        const PetscInt stencilSize = 13;
        /* The locations of the stencil points for the U equation */
        const PetscInt colU[] = {
          /*       */ i*M+Jp,
          im*M+Jp+1,  i*M+Jp+1,   ip*M+Jp+1,
          im*M+J,     i*M+J,      ip*M+J,
          im*M+J+1,               ip*M+J+1,
          /*       */ i*M+Jm,
          im*M+Jm+1,  i*M+Jm+1,   ip*M+Jm+1};
        /* The values at those points */
        PetscScalar valU[] = {
          /*               */ -c11/dy2,
          (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
          -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
          (c11-c10)/d4,                                       (c10-c11)/d4,
          /*               */ -c10/dy2,
          -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };

        /* The locations of the stencil points for the V equation */
        const PetscInt colV[] = {
          im*M+Jp,        i*M+Jp,     ip*M+Jp,
          /*           */ i*M+Jp+1,
          im*M+J,                     ip*M+J,
          im*M+J+1,       i*M+J+1,    ip*M+J+1,
          im*M+Jm,        i*M+Jm,     ip*M+Jm,
          /*           */ i*M+Jm+1 };
        /* The values at those points */
        PetscScalar valV[] = {
          (2*c11+c00)/d4,     (c00-c01)/d4,                   -(2*c11+c01)/d4,
          /*               */ -4*c11/dy2,
          2*(c11-c10)/d4,                                     2*(c10-c11)/d4,
          -c00/dx2,           4*(c11+c10)/dy2+(c01+c00)/dx2,  -c01/dx2,
          -(2*c10+c00)/d4,    (c01-c00)/d4,                   (2*c10+c01)/d4,
          /*               */ -4*c10/dy2 };
#endif
#endif

        /* Dragging ice experiences friction at the bed determined by the
        * basalDrag[x|y]() methods.  These may be a linear friction law or plastic. */
        if (intMask(mask[i][j]) == MASK_DRAGGING) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          valU[5] += basalDragx(beta, tauc, u, v, i, j);
          valV[7] += basalDragy(beta, tauc, u, v, i, j);
        }

        ierr = MatSetValues(A, 1, &rowU, stencilSize, colU, valU, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, stencilSize, colV, valV, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}


//! Computes the right-hand side of the linear problem for the SSA equations.
/*! 
@cond CONTINUUM
The right side of the SSA equations is just
   \f[ - \rho g H \nabla h \f]
because, in particular, the basal stress is put on the left side of the system.  
(See comment for assembleSSAMatrix().)  

@endcond
@cond NUMERIC
The surface slope is computed by centered difference onto the regular grid,
which may use periodic ghosting.  (Optionally the surface slope can be computed 
by only differencing into the grid for points at the edge; see test I.)

Note that the grid points with mask value MASK_SHEET correspond to the trivial 
equations
   \f[ \bar u_{ij} = \frac{uvbar[0][i-1][j] + uvbar[0][i][j]}{2}, \f]
and similarly for \f$\bar v_{ij}\f$.  That is, the vertically-averaged horizontal
velocity is already known for these points because it was computed (on the staggered
grid) using the SIA.
@endcond
 */
PetscErrorCode IceModel::assembleSSARhs(bool surfGradInward, Vec rhs) {
  // surfGradInward == true then differentiate h(x,y) inward from edge of grid,
  // so that certain solutions make sense on period grid; Test I, for now
  const PetscInt  Mx=grid.p->Mx, My=grid.p->My, M=2*My;
  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy;
  PetscErrorCode  ierr;
  PetscScalar     **mask, **h, **H, **uvbar[2];

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  /* matrix assembly loop */
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j;
      const PetscInt rowU = i*M + J;
      const PetscInt rowV = i*M + J+1;
      if (intMask(mask[i][j]) == MASK_SHEET) {
        ierr = VecSetValue(rhs, rowU, 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]),
                           INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(rhs, rowV, 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]),
                           INSERT_VALUES); CHKERRQ(ierr);
      } else {
        PetscScalar h_x, h_y;
        bool edge;
        if (surfGradInward) edge = ((i == 0) || (i == Mx-1) || (j == 0) || (j == My-1));
        if (surfGradInward && edge) {
          if (i == 0) {
            h_x = (h[i+1][j] - h[i][j]) / (dx);
            h_y = (h[i][j+1] - h[i][j-1]) / (2*dy);
          } else if (i == Mx-1) {
            h_x = (h[i][j] - h[i-1][j]) / (dx);
            h_y = (h[i][j+1] - h[i][j-1]) / (2*dy);
          } else if (j == 0) {
            h_x = (h[i+1][j] - h[i-1][j]) / (2*dx);
            h_y = (h[i][j+1] - h[i][j]) / (dy);
          } else if (j == My-1) {        
            h_x = (h[i+1][j] - h[i-1][j]) / (2*dx);
            h_y = (h[i][j] - h[i][j-1]) / (dy);
          } else {
            SETERRQ(1,"should not reach here: surfGradInward=TRUE & edge=TRUE but not at edge");
          }          
        } else { 
          h_x = (h[i+1][j] - h[i-1][j]) / (2*dx);
          h_y = (h[i][j+1] - h[i][j-1]) / (2*dy);
        }
        const PetscScalar r = ice.rho * grav * H[i][j];
        ierr = VecSetValue(rhs, rowU, -r*h_x, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(rhs, rowV, -r*h_y, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::moveVelocityToDAVectors(Vec x) {
  const PetscInt  M = 2 * grid.p->My;
  PetscErrorCode  ierr;
  PetscScalar     **u, **v, *uv;
  Vec             xLoc = SSAXLocal;

  /* Move the solution onto a grid which can be accessed normally. Since the parallel
   * layout of the vector x does not in general have anything to do with the DA based
   * vectors, we must scatter the entire vector to all processors. */
  ierr = VecScatterBegin(SSAScatterGlobalToLocal, x, xLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(SSAScatterGlobalToLocal, x, xLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGetArray(xLoc, &uv); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      u[i][j] = uv[i*M + 2*j];
      v[i][j] = uv[i*M + 2*j+1];
    }
  }
  ierr = VecRestoreArray(xLoc, &uv); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);

  // Communicate so that we have stencil width for geometry evolution
  ierr = DALocalToLocalBegin(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::updateNuViewers(Vec vNu[2], Vec vNuOld[2], bool updateNu_tView) {
  PetscErrorCode ierr;
  if (lognuView != PETSC_NULL) {
    PetscScalar  **nui, **nuj, **gg;  
    ierr = DAVecGetArray(grid.da2, g2, &gg); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vNu[0], &nui); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vNu[1], &nuj); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscReal avnu = 0.5 * (nui[i][j] + nuj[i][j]);
        if (avnu > 1.0e14) {
          gg[i][j] = log10(avnu);
        } else {
          gg[i][j] = 14.0;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vNu[0], &nui); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vNu[1], &nuj); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, g2, &gg); CHKERRQ(ierr);
    ierr = VecView(g2, lognuView); CHKERRQ(ierr);
  }
  if (nuView[0] != PETSC_NULL && nuView[1] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vNu[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, nuView[0]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vNu[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, nuView[1]); CHKERRQ(ierr);
  }
  if ((NuView[0] != PETSC_NULL) && (NuView[1] != PETSC_NULL) && updateNu_tView) {
    // note vNuOld[] contain *difference* of nu after testConvergenceofNu()
    ierr = DALocalToGlobal(grid.da2, vNuOld[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, NuView[0]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vNuOld[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, NuView[1]); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::velocitySSA() {
  PetscErrorCode ierr;

  Vec vNuDefault[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  ierr = velocitySSA(vNuDefault); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::velocitySSA(Vec vNu[2]) {
  // call this one directly if control over alloc of vNu[2] is needed (e.g. test J)
  PetscErrorCode ierr;
  KSP ksp = SSAKSP;
  Mat A = SSAStiffnessMatrix;
  Vec x = SSAX, rhs = SSARHS; // solve  A x = rhs
  Vec vNuOld[2] = {vWork2d[2], vWork2d[3]};
  Vec vubarOld = vWork2d[4], vvbarOld = vWork2d[5];
  PetscReal   norm, normChange, epsilon;
  PetscInt    its;
  KSPConvergedReason  reason;

  // We need to save the velocity from the last time step since we may have to
  // restart the iteration with larger values of epsilon.
  ierr = VecCopy(vubar, vubarOld); CHKERRQ(ierr);
  ierr = VecCopy(vvbar, vvbarOld); CHKERRQ(ierr);
  epsilon = ssaEpsilon;

  ierr = verbPrintf(4,grid.com, 
     "  [ssaEpsilon = %10.5e, ssaMaxIterations = %d\n",
     ssaEpsilon, ssaMaxIterations); CHKERRQ(ierr);
  ierr = verbPrintf(4,grid.com, 
     "   regularizingVelocitySchoof = %10.5e, regularizingLengthSchoof = %10.5e,\n",
     regularizingVelocitySchoof, regularizingLengthSchoof); CHKERRQ(ierr);
  ierr = verbPrintf(4,grid.com, 
     "   constantHardnessForSSA = %10.5e, ssaRelativeTolerance = %10.5e]\n",
    constantHardnessForSSA, ssaRelativeTolerance); CHKERRQ(ierr);

  for (PetscInt l=0; ; ++l) {
    ierr = computeEffectiveViscosity(vNu, epsilon); CHKERRQ(ierr);
    for (PetscInt k=0; k<ssaMaxIterations; ++k) {
      ierr = VecCopy(vNu[0], vNuOld[0]); CHKERRQ(ierr);
      ierr = VecCopy(vNu[1], vNuOld[1]); CHKERRQ(ierr);

      ierr = verbPrintf(3,grid.com, "  %d,%2d:", l, k); CHKERRQ(ierr);
      ierr = assembleSSAMatrix(vNu, A); CHKERRQ(ierr);
      ierr = assembleSSARhs((computeSurfGradInwardSSA == PETSC_TRUE), rhs); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com, "A:"); CHKERRQ(ierr);

      ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      /* -ksp_type will set it; defaults to GMRES(30)
       * -ksp_pc will set preconditioner; defaults to ILU */
      /* If you want to test different KSP methods, it may be helpful to
       * see how many iterations were necessary. Initial testing implies
       * that CGS takes roughly half the iterations of GMRES(30), but is 
       * not significantly faster. Furthermore, ILU and BJACOBI seem 
       * roughly equivalent. */
      //ierr = KSPSetType(ksp, KSPCGS); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
      ierr = KSPSolve(ksp, rhs, x); CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
      ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com, "S:%d,%d: ", its, reason); CHKERRQ(ierr);

      ierr = moveVelocityToDAVectors(x); CHKERRQ(ierr);
      ierr = computeEffectiveViscosity(vNu, epsilon); CHKERRQ(ierr);
      ierr = testConvergenceOfNu(vNu, vNuOld, &norm, &normChange); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com,"|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n",
                         norm, normChange/norm); CHKERRQ(ierr);

      /* extra diagnositic info */
      ierr = updateNuViewers(vNu, vNuOld, true); CHKERRQ(ierr);

      if (norm == 0 || normChange / norm < ssaRelativeTolerance) goto done;
    }
    if (epsilon > 0.0) {
       ierr = verbPrintf(1,grid.com,
                  "WARNING: Effective viscosity not converged after %d iterations\n"
                  "\twith epsilon=%8.2e. Retrying with epsilon * %8.2e.\n",
                  ssaMaxIterations, epsilon, DEFAULT_EPSILON_MULTIPLIER_SSA);
           CHKERRQ(ierr);
       ierr = VecCopy(vubarOld, vubar); CHKERRQ(ierr);
       ierr = VecCopy(vvbarOld, vvbar); CHKERRQ(ierr);
       epsilon *= DEFAULT_EPSILON_MULTIPLIER_SSA;
    } else {
       SETERRQ1(1, 
         "Effective viscosity not converged after %d iterations; epsilon=0.0.  Stopping.\n", 
         ssaMaxIterations);
    }
  }

  done:
  ierr = verbPrintf(3,grid.com, " "); CHKERRQ(ierr);  // has to do with summary appearance

  // check if user wants the matrix, rhs, and soln vector saved in Matlab readable format
  if (ssaSystemToASCIIMatlab == PETSC_TRUE) {
    ierr = writeSSAsystemMatlab(vNu); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeSSAsystemMatlab(Vec vNu[2]) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           file_name[PETSC_MAX_PATH_LEN], yearappend[PETSC_MAX_PATH_LEN];

  strcpy(file_name,ssaMatlabFilePrefix);
  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_%.0f.m", grid.p->year);
  strcat(file_name,yearappend);
  ierr = verbPrintf(2, grid.com, 
             "writing Matlab-readable file for SSA system A xsoln = rhs to file `%s' ...\n",
             file_name); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name);CHKERRQ(ierr);

  // save linear system; gives system A xsoln = rhs at last (nonlinear) iteration of SSA
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% PISM SSA linear system report.\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Writes matrix A (sparse) and vectors xsoln, rhs for the linear system\n"); 
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% (A xsoln = rhs) which was solved at the last step of the nonlinear iteration.\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Also writes the year, the coordinates x,y and their gridded versions xx,yy.\n"); 
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Also thickness (H), surface elevation (h), i-offset (staggered grid) \n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% vertically-integrated viscosity (nu_0 = nu * H), and j-offset version\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% of same thing (nu_1 = nu * H).\n\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% The thickness can be plotted by\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%%   >> imagesc(x,y,flipud(H')), axis square, colorbar \n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% for example. \n");
             CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\necho off\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAStiffnessMatrix,"A"); CHKERRQ(ierr);
  ierr = MatView(SSAStiffnessMatrix, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSARHS,"rhs"); CHKERRQ(ierr);
  ierr = VecView(SSARHS, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAX,"xsoln"); CHKERRQ(ierr);
  ierr = VecView(SSAX, viewer);CHKERRQ(ierr);

  // save coordinates (for viewing, primarily)
  ierr = PetscViewerASCIIPrintf(viewer,"year=%10.6f;\n",grid.p->year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "x=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",grid.p->Lx,grid.p->dx,grid.p->Lx);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "y=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",grid.p->Ly,grid.p->dy,grid.p->Ly);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n\n");  CHKERRQ(ierr);

  // also save thickness and effective viscosity
  ierr = VecViewDA2Matlab(vH, viewer, "H"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vh, viewer, "h"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vNu[0], viewer, "nu_0"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vNu[1], viewer, "nu_1"); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  return 0;
}


//! At all SSA points, update the velocity field.
PetscErrorCode IceModel::broadcastSSAVelocity() {
  // take ubar,vbar and update 3D horizontal velocities u,v
  // (w gets update later from vertVelocityFromIncompressibility())

  PetscErrorCode ierr;
  PetscScalar **mask, **ubar, **vbar, **ub, **vb;
  PetscScalar ***u, ***v, **uvbar[2];
  PetscScalar locCFLmaxdt2D = maxdt;
  
  /* This updates the 3D velocity field so that, for example, the temperature eqn
     knows about the velocity in the non-SHEET regions.  Basal vels also get updated. */
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  if (doSuperpose == PETSC_TRUE) {
    ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  }
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) != MASK_SHEET) {
        // apply glaciological-superposition-to-low-order if desired (and not floating)
        bool addVels = ((doSuperpose == PETSC_TRUE) && (modMask(mask[i][j]) == MASK_DRAGGING));
        
        for (PetscInt k=0; k<grid.p->Mz; ++k) {
          u[i][j][k] = (addVels) ? u[i][j][k] + ubar[i][j] : ubar[i][j];
          v[i][j][k] = (addVels) ? v[i][j][k] + vbar[i][j] : vbar[i][j];
        }
        ub[i][j] = (addVels) ? ub[i][j] + ubar[i][j] : ubar[i][j];
        vb[i][j] = (addVels) ? vb[i][j] + vbar[i][j] : vbar[i][j];
        
        if (addVels) { // now update ubar,vbar by adding SIA contribution
          ubar[i][j] += 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]);
          vbar[i][j] += 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]);
        }

        PetscScalar denom = PetscAbs(ubar[i][j])/grid.p->dx + PetscAbs(vbar[i][j])/grid.p->dy;
        denom += (0.01/secpera)/(grid.p->dx + grid.p->dy);  // make sure it's pos.
        locCFLmaxdt2D = PetscMin(locCFLmaxdt2D,1.0/denom);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);

  if (doSuperpose == PETSC_TRUE) {
    ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
    // now communicate modified ubar, vbar because modified
    ierr = DALocalToLocalBegin(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  }

  ierr = PetscGlobalMin(&locCFLmaxdt2D, &CFLmaxdt2D, grid.com); CHKERRQ(ierr);

  return 0;
}


//! At SSA points, correct the previously-computed basal frictional heating.
PetscErrorCode IceModel::correctBasalFrictionalHeating() {
  // recompute vRb in ice stream (MASK_DRAGGING) locations; zeros vRb in FLOATING
  PetscErrorCode  ierr;
  PetscScalar **ub, **vb, **mask, **Rb, **beta, **tauc;

  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        Rb[i][j] = 0.0;
      }
      if ((modMask(mask[i][j]) == MASK_DRAGGING) && (useSSAVelocity)) {
        const PetscScalar beta_x = basalDragx(beta, tauc, ub, vb, i, j);
        const PetscScalar beta_y = basalDragy(beta, tauc, ub, vb, i, j);
        const PetscScalar basal_stress_x = - beta_x * ub[i][j];
        const PetscScalar basal_stress_y = - beta_y * vb[i][j];
        // note next line uses *updated* ub,vb if doSuperpose == TRUE
        Rb[i][j] = - basal_stress_x * ub[i][j] - basal_stress_y * vb[i][j];
      } 
      // otherwise leave SIA-computed value alone
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  return 0;
}


//! At SSA points, correct the previously-computed volume strain-heating.
PetscErrorCode IceModel::correctSigma() {
  // recompute vSigma in ice stream and shelf (DRAGGING,FLOATING) locations
  PetscErrorCode  ierr;
  PetscScalar **H, **mask, ***Sigma, **ub, **vb, ***T;

  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);

  const PetscScalar dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  // next constant is the form of regularization used by C. Schoof 2006 "A variational
  // approach to ice streams" J Fluid Mech 556 pp 227--251
  const PetscReal  schoofReg = PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) != MASK_SHEET) {
        // hor. velocities don't depend on depth; use basal values
        const PetscScalar u_x = (ub[i+1][j] - ub[i-1][j])/(2*dx),
                          u_y = (ub[i][j+1] - ub[i][j-1])/(2*dy),
                          v_x = (vb[i+1][j] - vb[i-1][j])/(2*dx),
                          v_y = (vb[i][j+1] - vb[i][j-1])/(2*dy);
        const PetscScalar beta = PetscSqr(u_x) + PetscSqr(v_y)
                           + u_x * v_y + PetscSqr(0.5*(u_y + v_x));
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        const PetscScalar CC = 4 * beta / (ice.rho * ice.c_p);

        // apply glaciological-superposition-to-low-order if desired (and not floating)
        bool addVels = ((doSuperpose == PETSC_TRUE) && (modMask(mask[i][j]) == MASK_DRAGGING));

        for (PetscInt k=0; k<ks; ++k) {
          // use hydrostatic pressure; presumably this is not quite right in context 
          // of shelves and streams
          const PetscScalar pressure = ice.rho * grav * (H[i][j] - k * dz);
          const PetscScalar mvSigma = CC * ice.effectiveViscosity(schoofReg,
                                               u_x,u_y,v_x,v_y,T[i][j][k],pressure);
          Sigma[i][j][k] = (addVels) ? Sigma[i][j][k] + mvSigma : mvSigma;
        }
        for (PetscInt k=ks+1; k<grid.p->Mz; ++k) {
          Sigma[i][j][k] = 0.0;
        }
      }
      // otherwise leave SIA-computed value alone
    }
  }

  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  return 0;
}
