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

PetscErrorCode IceModel::velocity(bool updateVelocityAtDepth) {
  PetscErrorCode ierr;
  static PetscTruth firstTime = PETSC_TRUE;

  ierr = velocitySIAStaggered(!updateVelocityAtDepth); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  for (PetscInt k=4; k<8; ++k) { // communicate ub[o], vb[o] on staggered
    ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
  }

  ierr = basalSIAConditionsToRegular(); CHKERRQ(ierr);
  
  if (updateVelocityAtDepth) {

    for (PetscInt k=0; k<6; ++k) { // communicate I[o], J[o], Sigma[o] on staggered
      ierr = DALocalToLocalBegin(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
    }
    for (PetscInt k=0; k<4; ++k) { // communicate h_x[o], h_y[o] on staggered
      ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
    }
    for (PetscInt k=8; k<10; ++k) { // communicate Rb[o] on staggered
      ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
    }

    ierr = SigmaSIAToRegular(); CHKERRQ(ierr);
    ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
    ierr = verticalVelocitySIARegular(); CHKERRQ(ierr);

    if (noSpokesLevel > 0) {
      ierr = smoothSigma(); CHKERRQ(ierr); // communication occurs in here!
    }
  }

  if (useMacayealVelocity) { // communication happens within MacAyeal
    ierr = setupForMacayeal(DEFAULT_MINH_MACAYEAL,PETSC_FALSE); CHKERRQ(ierr);
    if (firstTime) {
      ierr = vertAveragedVelocityToRegular(); CHKERRQ(ierr);  // comm here
    }
    ierr = velocityMacayeal(); CHKERRQ(ierr); // comm here ...
    ierr = cleanupAfterMacayeal(DEFAULT_MINH_MACAYEAL); CHKERRQ(ierr);
    ierr = broadcastMacayealVelocity(); CHKERRQ(ierr);
    ierr = correctSigma(); CHKERRQ(ierr);
    ierr = correctBasalFrictionalHeating(); CHKERRQ(ierr);
  } else { // Note vertically averaged vels on regular grid (Vecs vubar and 
     // vvbar) are set by MacAyeal procedures above, including (in a reasonable 
     // manner) within the MASK_SHEET regions.  Also vubar and vvbar are used 
     // to initialize the next step of MacAyeal.  SO we want to NOT overwrite 
     // the just-computed regular grid velocities (as they will be needed to
     // initialize another step of MacAyeal; note first time exception).
     // BUT if we don't call MacAyeal in velocity() then
     // we need to get staggered onto regular for viewer and for summary().
    ierr = vertAveragedVelocityToRegular(); CHKERRQ(ierr); // comm here
  }

  if ((useMacayealVelocity) || (updateVelocityAtDepth)) {
    ierr = computeMaxVelocities(); CHKERRQ(ierr); // communication here for global max,min
  }

  firstTime = PETSC_FALSE;
  return 0;
}


PetscErrorCode IceModel::velocitySIAStaggered(bool faststep) {
  // Vertically-integrated velocities (i.e. vuvbar) and basal sliding velocities
  // (vub, vvb) are *always* updated.  Diffusivity is only updated by vertical 
  // integration if !faststep.
  // If (useIsothermalFlux) then vuvbar is computed by using
  // isothermalFlux_A_softness and isothermalFlux_n_exponent with the Glen
  // formula:
  // vuvbar = U = - 2 A (rho g)^n (n+2)^{-1} H^{n+1} |\grad h|^{n-1} \grad h;
  // Set (useIsothermalFlux) with setIsothermalFlux() or with `-isoflux'.

  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy, dz=grid.p->dz;
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub[2], **vb[2], **Rb[2];
  PetscScalar ***I[2], ***J[2], ***Sigma[2];
  PetscScalar **h, **H, **mask, **Df[2], **uvbar[2];
  PetscScalar ***T, ***gs;
  PetscScalar *delta, *K;
  const PetscScalar Gamma = (2.0 * isothermalFlux_A_softness
                             * pow(ice.rho * ice.grav, isothermalFlux_n_exponent)
                             / (isothermalFlux_n_exponent + 2.0));

  delta = new PetscScalar[grid.p->Mz];
  K = new PetscScalar[grid.p->Mz];

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vDf[0], &Df[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vDf[1], &Df[1]); CHKERRQ(ierr);

  // note basal temps get evaled by basal sliding law, even in faststep
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);

  if (!faststep) {
    ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[8], &Rb[0]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[9], &Rb[1]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
  }
  
  // staggered grid computation of: I, J, Sigma
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        // staggered point: o==0 is right, o==1 is up
        const PetscInt      oi = 1-o, oj=o;

        PetscScalar  slope;
        if (o==0) {     // If I-offset
          slope = (h[i+1][j] - h[i][j]) / dx;
          h_x[o][i][j] = slope;
          h_y[o][i][j] = (+ h[i+1][j+1] + h[i][j+1]
                          - h[i+1][j-1] - h[i][j-1]) / (4.0*dy);
        } else {        // J-offset
          slope = (h[i][j+1] - h[i][j]) / dy;
          h_y[o][i][j] = slope;
          h_x[o][i][j] = (+ h[i+1][j+1] + h[i+1][j]
                          - h[i-1][j+1] - h[i-1][j]) / (4.0*dx);
        }

        const PetscScalar   thickness = 0.5 * (H[i][j] + H[i+oi][j+oj]);
 
        if (thickness > 0) { 
          const PetscInt      ks = static_cast<PetscInt>(floor(thickness/dz));
          const PetscScalar   alpha =
                  sqrt(PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]));

          if (!faststep) { // will evaluate flow law at depth
            if (ks>grid.p->Mz) {
              ierr = PetscPrintf(grid.com,
                   "[[error LOCATION: i, j, ks, H = %5d %5d %5d %10.2f]]\n",i, j, ks, H[i][j]); 
              SETERRQ(1, "thickness overflow in SIA velocity: ks>Mz");
            }
            I[o][i][j][0] = 0; J[o][i][j][0] = 0; K[0] = 0;
            for (PetscInt k=0; k<=ks; ++k) {
              const PetscScalar   s = k * dz;
              const PetscScalar   pressure = ice.rho * ice.grav * (thickness-s);
              delta[k] = (2 * pressure * enhancementFactor
                          * ice.flow(alpha * pressure,
                                     0.5 * (T[i][j][k] + T[i+oi][j+oj][k]), pressure,
                                     0.5 * (gs[i][j][k] + gs[i+oi][j+oj][k])));
              // for Sigma, ignor mask value and assume SHEET; will be overwritten
              // by correctSigma() in iMmacayeal.cc
              Sigma[o][i][j][k] = delta[k] * PetscSqr(alpha) * pressure
                                    / (ice.rho * ice.c_p);
              if (k>0) { // trapezoid rule for I[][][][k] and K[k]
                I[o][i][j][k] = I[o][i][j][k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
                K[k] = K[k-1] + 0.5 * dz * ((s-dz)*delta[k-1] + s*delta[k]);
                J[o][i][j][k] = s * I[o][i][j][k] - K[k];
              }
            }
            for (PetscInt k=ks+1; k<grid.p->Mz; ++k) { // above the ice
              Sigma[o][i][j][k] = 0.0;
              I[o][i][j][k] = I[o][i][j][ks];
              J[o][i][j][k] = k * dz * I[o][i][j][ks]; // J[o][i][j][ks];
            }  

            // diffusivity for deformational flow (vs basal diffusivity, incorporated in ub,vb)
            if (useIsothermalFlux) {
              Df[o][i][j] = -Gamma * pow(thickness, isothermalFlux_n_exponent + 2)
                                        * pow(alpha, isothermalFlux_n_exponent - 1);
            } else { // usual thermocoupled case
              Df[o][i][j] = J[o][i][j][ks] + (thickness - ks*dz) * I[o][i][j][ks];
            }
            // basal frictional heating
            if (modMask(mask[i][j]) == MASK_FLOATING) {
              Rb[o][i][j] = 0.0;
            } else { // ignor ice streams; will be overwritten by
                     //   correctBasalFrictionalHeating() if useMacAyealVelocities==TRUE
              const PetscScalar P = ice.rho * ice.grav * thickness;
              const PetscScalar basal_stress_x = P * h_x[o][i][j];
              const PetscScalar basal_stress_y = P * h_y[o][i][j];
              Rb[o][i][j] = - basal_stress_x * ub[o][i][j] - basal_stress_y * vb[o][i][j];
            }
          } // done with evals, calcs at depth, and those that affect temp eqn

          // basal velocity
          const PetscScalar myx = -grid.p->Lx + grid.p->dx * i, 
                            myy = -grid.p->Ly + grid.p->dy * j,
                            myT = 0.5 * (T[i][j][0] + T[i+oi][j+oj][0]);
          const PetscScalar basalC =
            basal(myx, myy, thickness, myT, alpha, muSliding);
          ub[o][i][j] = - basalC * h_x[o][i][j];
          vb[o][i][j] = - basalC * h_y[o][i][j];
          // note (*slide) is either ub[o][i][j] or vb[o][i][j] as appropriate

          // vertically-averaged velocity; note uvbar[0][i][j] is  u  at right staggered
          // point (i+1/2,j) but uvbar[1][i][j] is  v  at up staggered point (i,j+1/2)
          // here we use stale (old) diffuivity if faststep
          PetscScalar  *slide;       
          if (o==0) {     // If I-offset
            slide = &ub[o][i][j];
          } else {        // J-offset
            slide = &vb[o][i][j];
          }
          uvbar[o][i][j] = - Df[o][i][j] * slope / thickness + (*slide);
          
        } else {  // zero thickness case
          ub[o][i][j] = 0;
          vb[o][i][j] = 0;
          uvbar[o][i][j] = 0;
          if (!faststep) {  // only zero these when would be updated anyway
            Df[o][i][j] = 0;
            Rb[o][i][j] = 0;
            for (PetscInt k=0; k < grid.p->Mz; k++) {
              Sigma[o][i][j][k] = 0;
              I[o][i][j][k] = 0;
              J[o][i][j][k] = 0;
            }
          }
        } 
      }
    }
  }  
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da2, vDf[0], &Df[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vDf[1], &Df[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

  if (!faststep) {
    ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[8], &Rb[0]); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[9], &Rb[1]); CHKERRQ(ierr);
  }
  delete [] delta; delete [] K;
  return 0;
}


PetscErrorCode IceModel::basalSIAConditionsToRegular() {
  // compute vub, vvb, vRb by averaging from staggered onto regular grid
  PetscErrorCode  ierr;
  PetscScalar **ubreg, **vbreg, **ub[2], **vb[2];

  ierr = DAVecGetArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ubreg); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vbreg); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // average basal vels onto regular grid
      ubreg[i][j] = 0.25 * (ub[0][i][j] + ub[0][i-1][j] +
                            ub[1][i][j] + ub[1][i][j-1]);
      vbreg[i][j] = 0.25 * (vb[0][i][j] + vb[0][i-1][j] +
                            vb[1][i][j] + vb[1][i][j-1]);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vub, &ubreg); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vbreg); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::SigmaSIAToRegular() {
  // average Sigma onto regular grid for use in the temperature equation
  PetscErrorCode  ierr;
  PetscScalar ***Sigma[2], **H, ***Sigmareg, **Rb[2], **Rbreg;

  ierr = DAVecGetArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigmareg); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[8], &Rb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[9], &Rb[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vRb, &Rbreg); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // average frictional heating onto regular grid
      Rbreg[i][j] = 0.25 * (Rb[0][i][j] + Rb[0][i-1][j] +
                            Rb[1][i][j] + Rb[1][i][j-1]);
      // horizontally average Sigma onto regular grid
      const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
      for (PetscInt k=0; k<ks; ++k) {
        Sigmareg[i][j][k] = 0.25 * (Sigma[0][i][j][k] + Sigma[0][i-1][j][k] +
                                  Sigma[1][i][j][k] + Sigma[1][i][j-1][k]);
      }
      for (PetscInt k=ks+1; k<grid.p->Mz; ++k) {
        Sigmareg[i][j][k] = 0.0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigmareg); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[8], &Rb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[9], &Rb[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vRb, &Rbreg); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::horizontalVelocitySIARegular() {
  // update regular grid horizontal velocities u,v at depth
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub, **vb;
  PetscScalar ***I[2], ***u, ***v;

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<grid.p->Mz; ++k) {
        u[i][j][k] =  ub[i][j] - 0.25 *
          (I[0][i][j][k]*h_x[0][i][j] + I[0][i-1][j][k]*h_x[0][i-1][j] +
           I[1][i][j][k]*h_x[1][i][j] + I[1][i][j-1][k]*h_x[1][i][j-1]);
        v[i][j][k] =  vb[i][j] - 0.25 *
          (I[0][i][j][k]*h_y[0][i][j] + I[0][i-1][j][k]*h_y[0][i-1][j] +
           I[1][i][j][k]*h_y[1][i][j] + I[1][i][j-1][k]*h_y[1][i][j-1]);
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::verticalVelocitySIARegular() {
  // update vertical velocity at depth
  
  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy, dz=grid.p->dz;
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub[2], **vb[2], **basalMeltRate;
  PetscScalar ***J[2], ***w;

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar baseDivU = (ub[0][i][j] - ub[0][i-1][j]) / dx
                                   + (vb[1][i][j] - vb[1][i][j-1]) / dy;
      for (PetscInt k=0; k<grid.p->Mz; ++k) {
        const PetscScalar wfromslide = - (k * dz) * baseDivU;
        const PetscScalar wmain =
                    (J[0][i][j][k]*h_x[0][i][j] - J[0][i-1][j][k]*h_x[0][i-1][j]) / dx
                  + (J[1][i][j][k]*h_y[1][i][j] - J[1][i][j-1][k]*h_y[1][i][j-1]) / dy;
        w[i][j][k] = wmain + wfromslide;
        if (includeBMRinContinuity == PETSC_TRUE) {
          w[i][j][k] -= capBasalMeltRate(basalMeltRate[i][j]);
        }
/* is the following really desirable? is w relative to base of ice or not?:
        if (modMask(mask[i][j]) != MASK_FLOATING) {
          // only add effect of sloped or moving bed if
          // ice is in contact with the bed
          const PetscScalar Ik = 0.25 *
            (I[0][i][j][k] + I[0][i-1][j][k] + I[1][i][j][k] + I[1][i][j-1][k]);
          const PetscScalar h_x0 = (h[i+1][j] - h[i-1][j]) / (2*dx);
          const PetscScalar h_y0 = (h[i][j+1] - h[i][j-1]) / (2*dy);
          const PetscScalar b_x0 = (b[i+1][j] - b[i-1][j]) / (2*dx);
          const PetscScalar b_y0 = (b[i][j+1] - b[i][j-1]) / (2*dy);
          w[i][j][k] += - (Ik * h_x0 - u[i][j][0]) * b_x0
                        - (Ik * h_y0 - v[i][j][0]) * b_y0;
          // DO NOT add dbed/dt here, because w is relative to bed
        }
*/
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceModel::capBasalMeltRate(const PetscScalar bMR) {
  const PetscScalar MAX_BASALMELTRATE = 0.5/secpera;
  if (bMR > MAX_BASALMELTRATE) {
    return MAX_BASALMELTRATE;
  } else if (bMR < -MAX_BASALMELTRATE) {
    return -MAX_BASALMELTRATE;
  } else {
    return bMR;
  }
}


PetscErrorCode IceModel::smoothSigma() {
  // does iterated smoothing of Sigma on regular grid; uses box stencil of width one
  PetscErrorCode  ierr;
  PetscScalar ***S, ***Snew, **H;
  // following is [0.35 0.30 0.35] outer product w itself; ref conversation with Orion L.:
  const PetscScalar c[3][3] = {{0.1225, 0.105,  0.1225},
                               {0.105,  0.09,   0.105},
                               {0.1225, 0.105,  0.1225}};

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);

  for (int count=0; count < noSpokesLevel; ++count) {
    ierr = DAVecGetArray(grid.da3, vSigma, &S); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[4], &Snew); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; ++k) {
          // note S[neighbor_i][neighbor_j][k] will be zero if outside of ice
          const PetscScalar SS =  c[0][0]*S[i-1][j-1][k] + c[0][1]*S[i-1][j][k] + c[0][2]*S[i-1][j+1][k]
            + c[1][0]*S[i][j-1][k]   + c[1][1]*S[i][j][k]   + c[1][2]*S[i][j+1][k]
            + c[2][0]*S[i+1][j-1][k] + c[2][1]*S[i+1][j][k] + c[2][2]*S[i+1][j+1][k];
          const PetscScalar myz = k * grid.p->dz;
          PetscScalar active = 0.0;  // build sum of coeffs for neighbors with ice at or above curr depth
          if (H[i-1][j-1] >= myz) { active += c[0][0]; } 
          if (H[i-1][j] >= myz) { active += c[0][1]; } 
          if (H[i-1][j+1] >= myz) { active += c[0][2]; } 
          if (H[i][j-1] >= myz) { active += c[1][0]; } 
          if (H[i][j+1] >= myz) { active += c[1][2]; } 
          if (H[i+1][j-1] >= myz) { active += c[2][0]; } 
          if (H[i+1][j] >= myz) { active += c[2][1]; } 
          if (H[i+1][j+1] >= myz) { active += c[2][2]; } 
          Snew[i][j][k] = SS / (active + c[1][1]);  // ensures all heating goes into ice not outside it
        }
        for (PetscInt k=ks+1; k<grid.p->Mz; ++k) {
          Snew[i][j][k] = 0;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da3, vSigma, &S); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Snew); CHKERRQ(ierr);
      
    // communicate ghosted values *and* transfer Snew to vSigma
    ierr = DALocalToLocalBegin(grid.da3, vWork3d[4], INSERT_VALUES, vSigma); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vWork3d[4], INSERT_VALUES, vSigma); CHKERRQ(ierr);
  } // for
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::vertAveragedVelocityToRegular() {
  // only 2D regular grid velocities are updated here
  
  PetscErrorCode ierr;
  PetscScalar **u, **v, **uvbar[2];

  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      u[i][j] = 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]);
      v[i][j] = 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]);
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::computeMaxVelocities() {
  // computes max velocities in 3D grid and also sets CFLmaxdt by CFL condition
  PetscErrorCode ierr;
  PetscScalar ***u, ***v, ***w, **H;
  PetscScalar locCFLmaxdt;

  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);

  locCFLmaxdt = maxdt;
  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt      ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      for (PetscInt k=0; k<ks; ++k) {
        const PetscScalar au = PetscAbs(u[i][j][k]);
        const PetscScalar av = PetscAbs(v[i][j][k]);
        maxu = PetscMax(maxu,au);
        maxv = PetscMax(maxv,av);
        PetscScalar tempdenom = (0.001/secpera)/(grid.p->dx + grid.p->dy);  // make sure it's pos.
        tempdenom += PetscAbs(au/grid.p->dx) + PetscAbs(av/grid.p->dy);
        if (!isMarginal) {
          const PetscScalar aw = PetscAbs(w[i][j][k]);
          maxw = PetscMax(maxw,aw);
          tempdenom += PetscAbs(aw/grid.p->dz);
        }
        locCFLmaxdt = PetscMin(locCFLmaxdt,1.0 / tempdenom); 
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxu, &gmaxu, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxv, &gmaxv, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxw, &gmaxw, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMin(&locCFLmaxdt, &CFLmaxdt, grid.com); CHKERRQ(ierr);

  return 0;
}
