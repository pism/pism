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

PetscErrorCode IceModel::velocity(bool updateSIAVelocityAtDepth) {
  PetscErrorCode ierr;
  static PetscTruth firstTime = PETSC_TRUE;

  ierr = velocitySIAstaggered(); CHKERRQ(ierr);
  
  if (updateSIAVelocityAtDepth) {
    ierr = velocitySIAregular(); CHKERRQ(ierr);
  }
  if (noSpokesLevel > 0) {
    ierr = smoothSigma(); CHKERRQ(ierr);
  }

  if (useMacayealVelocity) {
    ierr = setupForMacayeal(DEFAULT_MINH_MACAYEAL,PETSC_FALSE); CHKERRQ(ierr);
    if (firstTime) {
      ierr = mapStaggeredVelocityToStandard(); CHKERRQ(ierr);
    }
    ierr = velocityMacayeal(); CHKERRQ(ierr);
    ierr = cleanupAfterMacayeal(DEFAULT_MINH_MACAYEAL); CHKERRQ(ierr);
    ierr = broadcastMacayealVelocity(); CHKERRQ(ierr);
  } else { // note vertically averaged vels on standard grid (Vecs vubar and vvbar)
           // are set by MacAyeal procedures above,
           // including (in a reasonable manner) within the MASK_SHEET regions,
           // and we want to not overwrite the just-computed standard grid
           // velocities in order to initialize another step of MacAyeal,
           // BUT if we don't call the MacAyeal we need to get staggered onto 
           // standard for viewer and for summary()
    ierr = mapStaggeredVelocityToStandard(); CHKERRQ(ierr);
  }

  ierr = computeMaxVelocities(); CHKERRQ(ierr);

  firstTime = PETSC_FALSE;
  return 0;
}


PetscErrorCode IceModel::velocitySIAstaggered() {
  // Vertically-integrated velocities (i.e. vuvbar) are *always* updated.  If
  // (useIsothermalFlux) then vuvbar is computed by using
  // isothermalFlux_A_softness and isothermalFlux_n_exponent with the Glen
  // formula:
  // vuvbar = U = - 2 A (rho g)^n (n+2)^{-1} H^{n+1} |\grad h|^{n-1} \grad h;
  // Set (useIsothermalFlux) with setIsothermalFlux() or with `-isoflux'.

  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy, dz=grid.p->dz;
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub[2], **vb[2];
  PetscScalar ***I[2], ***J[2], ***Sigma[2];
  PetscScalar **h, **H, **mask, **uvbar[2];
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

  ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);

  /*
  * Offset grids, computation of diffusivity; computation of strain-heating term Sigma
  */
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        const PetscInt      oi = 1-o, oj=o;
        const PetscScalar   thickness = 0.5 * (H[i][j] + H[i+oi][j+oj]);
        PetscScalar         slope, *slide;
        if (o==0) {     // If I-offset
          slide = &ub[o][i][j];
          slope = (h[i+1][j] - h[i][j]) / dx;
          h_x[o][i][j] = slope;
          h_y[o][i][j] = (+ h[i+1][j+1] + h[i][j+1]
                          - h[i+1][j-1] - h[i][j-1]) / (4.0*dy);
        } else {        // J-offset
          slide = &vb[o][i][j];
          slope = (h[i][j+1] - h[i][j]) / dy;
          h_y[o][i][j] = slope;
          h_x[o][i][j] = (+ h[i+1][j+1] + h[i+1][j]
                          - h[i-1][j+1] - h[i-1][j]) / (4.0*dx);
        }

        if (thickness > 0) {
          const PetscScalar   alpha =
            sqrt(PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]));
          const PetscInt      ks = static_cast<PetscInt>(floor(thickness/dz));

          if (ks>grid.p->Mz) {
            SETERRQ(1, "Thickness overflow in SIA velocity: ks>Mz");
          }
          I[o][i][j][0] = 0; J[o][i][j][0] = 0; K[0] = 0;
          for (PetscInt k=0; k<=ks; ++k) {
            const PetscScalar   s = k * dz;
            const PetscScalar   pressure = ice.rho * ice.grav * (thickness-s);
            delta[k] = (2 * pressure * enhancementFactor
                        * ice.flow(alpha * pressure,
                                   0.5 * (T[i][j][k] + T[i+oi][j+oj][k]), pressure,
                                   0.5 * (gs[i][j][k] + gs[i+oi][j+oj][k])));
            Sigma[o][i][j][k] = (delta[k] * PetscSqr(alpha) * pressure
                                 / (ice.rho * ice.c_p));
            if (k>0) { // trapezoid rule for I[][][][k] and K[k]
              I[o][i][j][k] = I[o][i][j][k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
              K[k] = K[k-1] + 0.5 * dz * ((s-dz)*delta[k-1] + s*delta[k]);
              J[o][i][j][k] = s * I[o][i][j][k] - K[k];
            }
          }
          for (PetscInt k=ks+1; k<grid.p->Mz; ++k) { // Not all needed, but this is safe
            Sigma[o][i][j][k] = 0;
            I[o][i][j][k] = I[o][i][j][ks];
            J[o][i][j][k] = k * dz * I[o][i][j][ks]; // J[o][i][j][ks];
          }

          const PetscScalar myx = -grid.p->Lx + grid.p->dx * i, 
                            myy = -grid.p->Ly + grid.p->dy * j,
                            myT = 0.5 * (T[i][j][0] + T[i+oi][j+oj][0]);
          const PetscScalar basalC =
            basal(myx, myy, thickness, myT, alpha, muSliding);
          ub[o][i][j] = - basalC * h_x[o][i][j];
          vb[o][i][j] = - basalC * h_y[o][i][j];
          // note (*slide) is either ub[o][i][j] or vb[o][i][j] as appropriate
          if (useIsothermalFlux) {
            uvbar[o][i][j] = (-Gamma * pow(thickness, isothermalFlux_n_exponent + 1)
                              * pow(alpha, isothermalFlux_n_exponent - 1) * slope) +
              (*slide);
          } else { // usual thermocoupled case
            const PetscScalar D = J[o][i][j][ks] + (thickness - ks*dz) * I[o][i][j][ks];
            uvbar[o][i][j] = -D * slope / thickness + (*slide);
          }
        } else {  // zero thickness case:
          for (PetscInt k=0; k < grid.p->Mz; k++) {
            Sigma[o][i][j][k] = 0;
            I[o][i][j][k] = 0;
            J[o][i][j][k] = 0;
            ub[o][i][j] = 0;
            vb[o][i][j] = 0;
            uvbar[o][i][j] = 0;
          }
        } 
      }
    }
  }  
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

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

  ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  
  delete [] delta; delete [] K;
  return 0;
}


PetscErrorCode IceModel::velocitySIAregular() {
  // update velocities at depth
  
  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy, dz=grid.p->dz;
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub[2], **vb[2];
  PetscScalar ***I[2], ***J[2], ***Sigma[2];
  PetscScalar **H, **h, **b, **basalMeltRate, **mask;
  PetscScalar ***u, ***v, ***w, ***Sigma0;

  // next lines only make sense if called just after velocitySIAstaggered
  for (PetscInt k=0; k<6; ++k) {
    ierr = DALocalToLocalBegin(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
  }
  for (PetscInt k=0; k<8; ++k) {
    ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
  }

  // average Sigma onto a regular grid for use in the temperature equation;
  //   set Sigma to zero above ice
  ierr = DAVecGetArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) == MASK_SHEET) {
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; ++k) {
          Sigma0[i][j][k] = 0.25 * (Sigma[0][i][j][k] + Sigma[0][i-1][j][k] +
                                    Sigma[1][i][j][k] + Sigma[1][i][j-1][k]);
        }
        for (PetscInt k=ks+1; k<grid.p->Mz; ++k) {
          Sigma0[i][j][k] = 0.0;
        }
      } else { // add ocean heat flux to bottom layer Sigma on ice *shelves*, but 
        // otherwise set Sigma to zero on ice shelves and ice streams (MacAyeal)
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          Sigma0[i][j][0] = DEFAULT_OCEAN_HEAT_FLUX / (ice.rho * ice.c_p * grid.p->dz);
        } else {
          Sigma0[i][j][0] = 0.0;  // no heating in streams at all
        }
        for (PetscInt k=1; k<grid.p->Mz; ++k) {
          Sigma0[i][j][k] = 0.0;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma0); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  // Compute velocity at depth
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // compute u and v at depth
      const PetscScalar   u0 = 0.25 *
        (ub[0][i][j] + ub[0][i-1][j] + ub[1][i][j] + ub[1][i][j-1]);
      const PetscScalar   v0 = 0.25 *
        (vb[0][i][j] + vb[0][i-1][j] + vb[1][i][j] + vb[1][i][j-1]);
      for (PetscInt k=0; k<grid.p->Mz; ++k) {
        u[i][j][k] =  u0 - 0.25 *
          (I[0][i][j][k]*h_x[0][i][j] + I[0][i-1][j][k]*h_x[0][i-1][j] +
           I[1][i][j][k]*h_x[1][i][j] + I[1][i][j-1][k]*h_x[1][i][j-1]);
        v[i][j][k] =  v0 - 0.25 *
          (I[0][i][j][k]*h_y[0][i][j] + I[0][i-1][j][k]*h_y[0][i-1][j] +
           I[1][i][j][k]*h_y[1][i][j] + I[1][i][j-1][k]*h_y[1][i][j-1]);
      }

      // vertical velocity
      const PetscScalar wmelt = - basalMeltRate[i][j];  // contributed by temperature model
      for (PetscInt k=0; k<grid.p->Mz; ++k) {
        const PetscScalar wfromslide = - (k * dz)
                    * (+ (ub[0][i][j] - ub[0][i-1][j]) / dx
                       + (vb[1][i][j] - vb[1][i][j-1]) / dy);
        const PetscScalar wmain =
                    (J[0][i][j][k]*h_x[0][i][j] - J[0][i-1][j][k]*h_x[0][i-1][j]) / dx
                  + (J[1][i][j][k]*h_y[1][i][j] - J[1][i][j-1][k]*h_y[1][i][j-1]) / dy;
        w[i][j][k] = wmain + wfromslide + wmelt;
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
          // FIXME: add dbed/dt here?  or not, because w is relative to bed
        }        
      } // for k

    }
  }

  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);

  return 0;
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


PetscErrorCode IceModel::mapStaggeredVelocityToStandard() {
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

  locCFLmaxdt = maxdt * tempskip;
  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt      ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
      for (PetscInt k=0; k<ks; ++k) {
        const PetscScalar au = PetscAbs(u[i][j][k]);
        const PetscScalar av = PetscAbs(v[i][j][k]);
        const PetscScalar aw = PetscAbs(w[i][j][k]);
        maxu = PetscMax(maxu,au);
        maxv = PetscMax(maxv,av);
        maxw = PetscMax(maxw,aw);
        const PetscScalar tempdenom = PetscAbs(au/grid.p->dx) + PetscAbs(av/grid.p->dy)
          + PetscAbs(aw/grid.p->dz);
        locCFLmaxdt = PetscMin(locCFLmaxdt,1.0 / tempdenom); 
      }
    }
  }
  ierr = PetscGlobalMax(&maxu, &gmaxu, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxv, &gmaxv, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxw, &gmaxw, grid.com); CHKERRQ(ierr);
  
  ierr = PetscGlobalMin(&locCFLmaxdt, &CFLmaxdt, grid.com); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  return 0;
}
