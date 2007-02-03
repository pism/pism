// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"


PetscErrorCode IceModel::setupForMacayeal(const PetscScalar minH, const PetscTruth adjustMask) {
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
      // FIXME: at this point we should probably not call with adjustMask == 
      // PETSC_TRUE because this cutoff is (I think) kludgier than the vote-by-
      // neighbors scheme in IceModel::updateSurfaceElevationandMask() and the
      // initialization of the mask using the mass balance velocity values (see
      // iMIOnetcdf.cc).
      if (adjustMask==PETSC_TRUE) {
        /* When this works properly, we should be able to tune it to reproduce ice
        * streams. It will only function correctly after the temperature has come to
        * equilibrium.                                                          */
        const PetscScalar speedSIA = sqrt(PetscSqr(u[i][j][0]) + PetscSqr(v[i][j][0]));
        const PetscScalar speed = sqrt(PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]));
        if (intMask(mask[i][j]) == MASK_SHEET && speedSIA > DEFAULT_MIN_SHEET_TO_DRAGGING/secpera) {
          mask[i][j] = MASK_DRAGGING;
        } else if (intMask(mask[i][j]) == MASK_DRAGGING
                   && speed < DEFAULT_MAX_SPEED_DRAGGING_TO_SHEET/secpera
                   && speedSIA < DEFAULT_MAX_SPEEDSIA_DRAGGING_TO_SHEET/secpera) {
          mask[i][j] = MASK_SHEET;
        }
      }
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
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::cleanupAfterMacayeal(const PetscScalar minH) {
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

PetscErrorCode IceModel::computeEffectiveViscosity(Vec vNu[2], PetscReal epsilon) {
  PetscErrorCode ierr;

  // 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of 30 MPa yr
  if (useConstantNuForMacAyeal == PETSC_TRUE) {
    ierr = VecSet(vNu[0], constantNuForMacAyeal); CHKERRQ(ierr);
    ierr = VecSet(vNu[1], constantNuForMacAyeal); CHKERRQ(ierr);
    return 0;
  } else { // initialize for below
    ierr = VecSet(vNu[0], 0.0); CHKERRQ(ierr);
    ierr = VecSet(vNu[1], 0.0); CHKERRQ(ierr);
  }
  
  /*
  * We need to compute integrated effective viscosity. It is locally determined by the
  * strain rates and temperature field.
  */
  PetscScalar **mask, **H, ***T, **nu[2], **u, **v;
  
  // next constant is the form of regularization used by C. Schoof 2006 "A variational
  // approach to ice streams" J Fluid Mech 556 pp 227--251
  const PetscReal  schoofReg = PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof);
  
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
          if (useConstantHardnessForMacAyeal == PETSC_FALSE) { // usual temperature-dependent case
            // "nu" is really "nu H"!
            nu[o][i][j] = ice.effectiveViscosityColumn(schoofReg,
                                    myH,dz, u_x, u_y, v_x, v_y, T[i][j], T[i+oi][j+oj]);
          } else { // constant \bar B case, i.e for EISMINT ROSS
            // "nu" is really "nu H"!
            nu[o][i][j] = myH * constantHardnessForMacAyeal * 0.5 *
              pow(schoofReg + PetscSqr(u_x) + PetscSqr(v_y) + 0.25*PetscSqr(u_y+v_x) + u_x*v_y,
                  -(1.0/3.0));
          }

          if (! finite(nu[o][i][j]) || false) {
            ierr = PetscPrintf(grid.com, "nu[%d][%d][%d] = %e\n", o, i, j, nu[o][i][j]);
            CHKERRQ(ierr); 
            ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n", u_x, u_y, v_x, v_y);
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

/* use -constant_nu 300.0 if you want about these values:
  ierr = VecSet(vNu[0], 1.0e16); CHKERRQ(ierr);
  ierr = VecSet(vNu[1], 1.0e16); CHKERRQ(ierr);
*/
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


PetscErrorCode IceModel::assembleMacayealMatrix(Vec vNu[2], Mat A) {
  const PetscInt  Mx=grid.p->Mx, My=grid.p->My, M=2*My;
  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy;
  const PetscScalar   one = 1.0;
  PetscErrorCode  ierr;
  PetscScalar     **mask, **nu[2], **u, **v;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  /* matrix assembly loop */
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
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
        * basalDrag() method. This may be a linear friction law. */
        if (intMask(mask[i][j]) == MASK_DRAGGING) {
          // We do dragging implicitly now.
          // use the value from Hulbe & MacAyeal (1999), p. 25,356
          valU[5] += basalDragx(u, v, i, j);
          valV[7] += basalDragy(u, v, i, j);
        }

        ierr = MatSetValues(A, 1, &rowU, stencilSize, colU, valU, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, stencilSize, colV, valV, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::assembleMacayealRhs(bool surfGradInward, Vec rhs) {
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
#if 0
        /* This equation seems to behave badly in the event of large surface
        * slopes. To handle this, we cap the slope for the purpose of this
        * computation only. */
        const PetscScalar   maxSlope = DEFAULT_MAXSLOPE_MACAYEAL;
        if (PetscAbsScalar(h_x) > maxSlope) h_x = h_x * maxSlope / PetscAbsScalar(h_x);
        if (PetscAbsScalar(h_y) > maxSlope) h_y = h_y * maxSlope / PetscAbsScalar(h_y);
#endif
        const PetscScalar r = ice.rho * ice.grav * H[i][j];  // CHANGE OF SIGN
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
  Vec             xLoc = MacayealXLocal;

  /* Move the solution onto a grid which can be accessed normally. Since the parallel
  * layout of the vector x does not in general have anything to do with the DA based
  * vectors, we must scatter the entire vector to all processors. */
  ierr = VecScatterBegin(x, xLoc, INSERT_VALUES, SCATTER_FORWARD,
                         MacayealScatterGlobalToLocal); CHKERRQ(ierr);
  ierr = VecScatterEnd(x, xLoc, INSERT_VALUES, SCATTER_FORWARD,
                       MacayealScatterGlobalToLocal); CHKERRQ(ierr);
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

PetscErrorCode IceModel::velocityMacayeal() {
//  PetscInt maxIterations;
  PetscErrorCode ierr;
  KSP ksp = MacayealKSP;
  Mat A = MacayealStiffnessMatrix;
  Vec x = MacayealX, rhs = MacayealRHS; // solve  A x = rhs
  Vec vNu[2] = {vWork2d[0], vWork2d[1]};
  Vec vNuOld[2] = {vWork2d[2], vWork2d[3]};
  Vec vubarOld = vWork2d[4], vvbarOld = vWork2d[5];
  PetscReal   norm, normChange, epsilon;
  PetscInt    its;
  KSPConvergedReason  reason;

  // We need to save the velocity from the last time step since we may have to
  // restart the iteration with larger values of epsilon.
  ierr = VecCopy(vubar, vubarOld); CHKERRQ(ierr);
  ierr = VecCopy(vvbar, vvbarOld); CHKERRQ(ierr);
  epsilon = macayealEpsilon;
//  ierr = verbPrintf(3,grid.com, "  iteration to solve ice stream/shelf equations:\n"); CHKERRQ(ierr);

  ierr = verbPrintf(5,grid.com, 
     "  [macayealEpsilon = %10.5e, EPSILON_MULTIPLIER_MACAYEAL = %6.2e, macayealMaxIterations = %d\n",
     macayealEpsilon, DEFAULT_EPSILON_MULTIPLIER_MACAYEAL,macayealMaxIterations); CHKERRQ(ierr);
  ierr = verbPrintf(5,grid.com, 
     "   regularizingVelocitySchoof = %10.5e, regularizingLengthSchoof = %10.5e,\n",
     regularizingVelocitySchoof, regularizingLengthSchoof); CHKERRQ(ierr);
  ierr = verbPrintf(5,grid.com, 
     "   constantHardnessForMacAyeal = %10.5e, macayealRelativeTolerance = %10.5e]\n",
    constantHardnessForMacAyeal, macayealRelativeTolerance); CHKERRQ(ierr);

  for (PetscInt l=0; ; ++l) {
    ierr = computeEffectiveViscosity(vNu, epsilon); CHKERRQ(ierr);
    for (PetscInt k=0; k<macayealMaxIterations; ++k) {
      ierr = VecCopy(vNu[0], vNuOld[0]); CHKERRQ(ierr);
      ierr = VecCopy(vNu[1], vNuOld[1]); CHKERRQ(ierr);

      ierr = verbPrintf(3,grid.com, "  %d,%2d:", l, k); CHKERRQ(ierr);
      ierr = assembleMacayealMatrix(vNu, A); CHKERRQ(ierr);
      ierr = assembleMacayealRhs((computeSurfGradInwardMacAyeal == PETSC_TRUE), rhs); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com, "A:"); CHKERRQ(ierr);

#if 0
        PetscReal vec_norm, mat_norm;
        ierr = MatNorm(A, NORM_FROBENIUS, &mat_norm); CHKERRQ(ierr);
        ierr = verbPrintf(3,grid.com, "\nNorm(A) = %e\n", mat_norm); CHKERRQ(ierr);
        ierr = MatMult(A, rhs, x); CHKERRQ(ierr);
        ierr = VecNorm(x, NORM_2, &vec_norm); CHKERRQ(ierr);
        ierr = verbPrintf(3,grid.com, "Norm(A * rhs) = % e\n", vec_norm); CHKERRQ(ierr);
#endif

      ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      /* -ksp_type will set it; defaults to GMRES(30)
       * -ksp_pc will set preconditioner; defaults to ILU
       * If you want to test different KSP methods, it may be helpful to
       * see how many iterations were necessary. Initial testing implies
       * that CGS takes roughly half the iterations of GMRES(30)
       * [default], but is not significantly faster. Furthermore, ILU
       * [default] and BJACOBI seem roughly equivalent. */
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

      if (norm == 0 || normChange / norm < macayealRelativeTolerance) goto done;
    }
    if (epsilon > 0.0) {
       ierr = verbPrintf(1,grid.com,
                  "WARNING: Effective viscosity not converged after %d iterations\n"
                  "\twith epsilon=%8.2e. Retrying with epsilon * %8.2e.\n",
                  macayealMaxIterations, epsilon, DEFAULT_EPSILON_MULTIPLIER_MACAYEAL);
           CHKERRQ(ierr);
       ierr = VecCopy(vubarOld, vubar); CHKERRQ(ierr);
       ierr = VecCopy(vvbarOld, vvbar); CHKERRQ(ierr);
       epsilon *= DEFAULT_EPSILON_MULTIPLIER_MACAYEAL;
    } else {
       SETERRQ1(1, "Effective viscosity not converged after %d iterations; epsilon=0.0.  Stopping.\n", 
                macayealMaxIterations);
    }
  }

  done:
  ierr = verbPrintf(3,grid.com, " "); CHKERRQ(ierr);  // has to do with summary appearance
  return 0;
}


PetscErrorCode IceModel::broadcastMacayealVelocity() {
  PetscErrorCode ierr;
  PetscScalar **mask, **b, **basalMeltRate, **ubar, **vbar, **ub, **vb;
  PetscScalar ***u, ***v, ***w;

  /* This allows the temperature equation to know about the velocity in the non-SHEET
  * regions. It is unnecessary if the temperature equation pays attention to the mask,
  * thus upwinding the 2-D velocity fields for the advection terms in the non-SHEET
  * regions.  Note that basal velocities also get updated. */
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) != MASK_SHEET) {
        ub[i][j] = ubar[i][j];
        vb[i][j] = vbar[i][j];
        for (PetscInt k=0; k<grid.p->Mz; ++k) {
          u[i][j][k] = ubar[i][j];
          v[i][j][k] = vbar[i][j];

#if 0
          // Upwind the horizontal velocities for the purpose of calculating
          // velocity gradients.
          const PetscScalar u_x = ((ubar[i][j] < 0)
                                   ? (ubar[i+1][j] - ubar[i][j]) / grid.p->dx
                                   : (ubar[i][j] - ubar[i-1][j]) / grid.p->dx);
          const PetscScalar v_y = ((vbar[i][j] < 0)
                                   ? (vbar[i][j+1] - vbar[i][j]) / grid.p->dy
                                   : (vbar[i][j] - vbar[i][j-1]) / grid.p->dy);
#else
          // ELB says: I see no reason to upwind in this context
          const PetscScalar u_x = (ubar[i+1][j] - ubar[i-1][j]) / (2.0*grid.p->dx);
          const PetscScalar v_y = (vbar[i][j+1] - vbar[i][j-1]) / (2.0*grid.p->dy);
#endif
          w[i][j][k] = - k * grid.p->dz * (u_x + v_y) 
                       - basalMeltRate[i][j];
          
          if (intMask(mask[i][j]) == MASK_DRAGGING) { // add velocity contribution from
                                 // contact with sloped and/or moving bed
            const PetscScalar   b_x0 = (b[i+1][j] - b[i-1][j]) / (2 * grid.p->dx);
            const PetscScalar   b_y0 = (b[i][j+1] - b[i][j-1]) / (2 * grid.p->dy);
            w[i][j][k] +=  + u[i][j][k] * b_x0 + v[i][j][k] * b_y0
                           + 0;  // db/dt  FIXME?: use uplift? probably NO:
                                 //   vert velocity is relative!
          }
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::correctBasalFrictionalHeating() {
  // recompute vRb in ice stream (MASK_DRAGGING) locations; zeros vRb in FLOATING
  PetscErrorCode  ierr;
  PetscScalar **ub, **vb, **mask, **Rb;

  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        Rb[i][j] = 0.0;
      }
      if ((modMask(mask[i][j]) == MASK_DRAGGING) && (useMacayealVelocity)) {
        const PetscScalar beta = basalDrag(ub[i][j], vb[i][j]);
        const PetscScalar basal_stress_x = - beta * ub[i][j];
        const PetscScalar basal_stress_y = - beta * vb[i][j];
        Rb[i][j] = - basal_stress_x * ub[i][j] - basal_stress_y * vb[i][j];
      } 
      // otherwise leave SIA-computed value alone
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  return 0;
}


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
      int m = modMask(mask[i][j]);
      if (( (m == MASK_DRAGGING) || (m == MASK_FLOATING) ) && (useMacayealVelocity)) {
        // hor. velocities don't depend on depth; use basal values
        const PetscScalar u_x = (ub[i+1][j] - ub[i-1][j])/(2*dx),
                          u_y = (ub[i][j+1] - ub[i][j-1])/(2*dy),
                          v_x = (vb[i+1][j] - vb[i-1][j])/(2*dx),
                          v_y = (vb[i][j+1] - vb[i][j-1])/(2*dy);
        const PetscScalar beta = PetscSqr(u_x) + PetscSqr(v_y)
                           + u_x * v_y + PetscSqr(0.5*(u_y + v_x));
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        const PetscScalar CC = 4 * beta / (ice.rho * ice.c_p);
        for (PetscInt k=0; k<ks; ++k) {
          // use hydrostatic pressure; presumably this is not quite right in context 
          // of shelves and streams
          const PetscScalar pressure = ice.rho * ice.grav * (H[i][j] - k * dz);
          Sigma[i][j][k] = CC * ice.effectiveViscosity(schoofReg,
                                         u_x,u_y,v_x,v_y,T[i][j][k],pressure);;
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
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  return 0;
}
