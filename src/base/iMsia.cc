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
#include <petscda.h>
#include "iceModel.hh"


//! Compute the surface gradient in advance of the SIA velocity computation.
/*! 
There are two methods for computing the surface gradient.  The default is to transform the 
thickness to something more regular and differentiate that.  In particular, as shown 
in (Calvo et al 2002) for the flat bed and \f$n=3\f$ case, if we define
	\f[\eta = H^{(2n+2)/n}\f]
then \f$\eta\f$ is more regular near the margin than \f$H\f$.  So the default method for computing
the surface gradient is to compute
   \f[\nabla h = \frac{n}{(2n+2)} \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b,\f]
recalling that \f$h = H + b\f$.  This method is only applied when \f$\eta > 0\f$ at a given point;
otherwise \f$\nabla h = \nabla b\f$.

We are computing this gradient by finite differences onto a staggered grid.  We do so 
by centered differences using (roughly) the same method for \f$\eta\f$ and 
\f$b\f$ that (Mahaffy 1976) applies directly to the surface elevation \f$h\f$.

The optional method is to directly differentiate the surface elevation \f$h\f$ by the
(Mahaffy 1976) method.
 */
PetscErrorCode IceModel::surfaceGradientSIA() {
  PetscErrorCode  ierr;

  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy;
  PetscScalar **h_x[2], **h_y[2];

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);

  if (transformForSurfaceGradient == PETSC_TRUE) {
    PetscScalar **eta, **b, **H;
    const PetscScalar n = isothermalFlux_n_exponent; // presumably 3.0
    const PetscScalar etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                      invpow  = 1.0 / etapow,
                      dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
    // compute eta = H^{8/3}, which is more regular, on reg grid
    ierr = DAVecGetArray(grid.da2, vWork2d[8], &eta); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        eta[i][j] = pow(H[i][j], etapow);
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[8], &eta); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    // communicate eta: other processors will need ghosted for d/dx and d/dy
    ierr = DALocalToLocalBegin(grid.da2, vWork2d[8], INSERT_VALUES, vWork2d[8]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vWork2d[8], INSERT_VALUES, vWork2d[8]); CHKERRQ(ierr);
    // now use Mahaffy on eta to get grad h on staggered;
    // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b
    ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[8], &eta); CHKERRQ(ierr);
    for (PetscInt o=0; o<2; o++) {
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
          if (o==0) {     // If I-offset
            const PetscScalar mean_eta = 0.5 * (eta[i+1][j] + eta[i][j]);
            if (mean_eta > 0.0) {
              const PetscScalar factor = invpow * pow(mean_eta, dinvpow);
              h_x[o][i][j] = factor * (eta[i+1][j] - eta[i][j]) / dx;
              h_y[o][i][j] = factor * (+ eta[i+1][j+1] + eta[i][j+1]
                                     - eta[i+1][j-1] - eta[i][j-1]) / (4.0*dy);
            } else {
              h_x[o][i][j] = 0.0;
              h_y[o][i][j] = 0.0;
            }
            // now add bed slope to get actual h_x,h_y
            h_x[o][i][j] += (b[i+1][j] - b[i][j]) / dx;
            h_y[o][i][j] += (+ b[i+1][j+1] + b[i][j+1]
                             - b[i+1][j-1] - b[i][j-1]) / (4.0*dy);
          } else {        // J-offset
            const PetscScalar mean_eta = 0.5 * (eta[i][j+1] + eta[i][j]);
            if (mean_eta > 0.0) {
              const PetscScalar factor = invpow * pow(mean_eta, dinvpow);
              h_y[o][i][j] = factor * (eta[i][j+1] - eta[i][j]) / dy;
              h_x[o][i][j] = factor * (+ eta[i+1][j+1] + eta[i+1][j]
                                     - eta[i-1][j+1] - eta[i-1][j]) / (4.0*dx);
            } else {
              h_y[o][i][j] = 0.0;
              h_x[o][i][j] = 0.0;
            }
            // now add bed slope to get actual h_x,h_y
            h_y[o][i][j] += (b[i][j+1] - b[i][j]) / dy;
            h_x[o][i][j] += (+ b[i+1][j+1] + b[i+1][j]
                             - b[i-1][j+1] - b[i-1][j]) / (4.0*dx);
          }
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[8], &eta); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  } else {  // if !transformForSurfaceGradient; the old way
    PetscScalar **h;
    ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
    for (PetscInt o=0; o<2; o++) {
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
          if (o==0) {     // If I-offset
            h_x[o][i][j] = (h[i+1][j] - h[i][j]) / dx;
            h_y[o][i][j] = (+ h[i+1][j+1] + h[i][j+1]
                            - h[i+1][j-1] - h[i][j-1]) / (4.0*dy);
          } else {        // J-offset
            h_y[o][i][j] = (h[i][j+1] - h[i][j]) / dy;
            h_x[o][i][j] = (+ h[i+1][j+1] + h[i+1][j]
                            - h[i-1][j+1] - h[i-1][j]) / (4.0*dx);
          }
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  } // end if (transformForSurfaceGradient)
  
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  return 0;
}


//!  Compute the vertically-averaged horizontal velocity according to the SIA.
/*!
See the comment for massBalExplicitStep() before reading the rest of this comment.

Note that one may write 
  \f[ \mathbf{q} = \bar{\mathbf{U}} H = D \nabla h + \mathbf{U}_b \cdot H\f]
in shallow ice approximation (SIA) areas.  Here \f$h\f$ is the surface elevation of the ice
\f$\mathbf{U}_b\f$ is the basal sliding velocity, and \f$D\f$ is the diffusivity (which 
is computed in this method).  At the end of this routine the value of the vertically-averaged
horizontal velocity \f$\bar{\mathbf{U}}\f$ is known at all staggered grid points.

The surface slope \f$\nabla h\f$ is needed on the staggered grid although the surface 
elevation \f$h\f$ itself is known on the regular grid.  The scheme used for this is
the one first proposed in the context of ice sheets by Mahaffy (1976).  That is, the method 
is "type I" in the classification described in (Hindmarsh and Payne 1996).

This routine also calls the part of the basal dynamical model applicable to the SIA; see 
basalVelocity().  The basal sliding velocity is computed for all SIA 
points.  This routine also computes the basal frictional heating and the volume
strain-heating.  Note that the SIA is used at all points on the grid in this routine 
but that the resulting vertically-averaged horizontal velocity is overwritten by different 
values at SSA points.  See correctBasalFrictionalHeating() and correctSigma().
 */
PetscErrorCode IceModel::velocitySIAStaggered(bool faststep) {
  // Vertically-integrated velocities (i.e. vuvbar) and basal sliding velocities
  // (vub, vvb) are *always* updated.  Diffusivity is only updated by vertical 
  // integration if !faststep.
  // If (useIsothermalFlux) then vuvbar is computed by using
  // isothermalFlux_A_softness and isothermalFlux_n_exponent with the Glen
  // formula:
  // vuvbar = U = - 2 A (rho g)^n (n+2)^{-1} H^{n+1} |\grad h|^{n-1} \grad h;
  // Set (useIsothermalFlux) with setIsothermalFlux() or with `-isoflux'.

  const PetscScalar   dz=grid.p->dz;
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub[2], **vb[2];
  PetscScalar ***I[2], ***J[2], ***Sigma[2];
  PetscScalar **h, **H, **Df[2], **uvbar[2];
  PetscScalar ***T, ***gs;
  PetscScalar *delta, *K;
  const PetscScalar Gamma = (2.0 * isothermalFlux_A_softness
                             * pow(ice.rho * grav, isothermalFlux_n_exponent)
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
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vDf[0], &Df[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vDf[1], &Df[1]); CHKERRQ(ierr);

  // note basal temps and melt rate get evaled by basal sliding law, even in faststep
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
  
  // staggered grid computation of: I, J, Sigma
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        // staggered point: o==0 is right, o==1 is up
        const PetscInt      oi = 1-o, oj=o;

        PetscScalar  slope;
        if (o==0) {     // If I-offset
          slope = h_x[o][i][j];
        } else {        // J-offset
          slope = h_y[o][i][j];
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
              const PetscScalar   pressure = ice.rho * grav * (thickness-s);
              delta[k] = (2 * pressure * enhancementFactor
                          * ice.flow(alpha * pressure,
                                     0.5 * (T[i][j][k] + T[i+oi][j+oj][k]), pressure,
                                     0.5 * (gs[i][j][k] + gs[i+oi][j+oj][k])));
              // for Sigma, ignor mask value and assume SHEET; will be overwritten
              // by correctSigma() in iMssa.cc
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
          } // done with evals, calcs at depth, and those that affect temp eqn

          // basal velocity
          const PetscScalar myx = -grid.p->Lx + grid.p->dx * i, 
                            myy = -grid.p->Ly + grid.p->dy * j,
                            myT = 0.5 * (T[i][j][0] + T[i+oi][j+oj][0]);
          const PetscScalar basalC =
            basalVelocity(myx, myy, thickness, myT, alpha, muSliding);
          ub[o][i][j] = - basalC * h_x[o][i][j];
          vb[o][i][j] = - basalC * h_y[o][i][j];

          // vertically-averaged velocity; note uvbar[0][i][j] is  u  at right staggered
          // point (i+1/2,j) but uvbar[1][i][j] is  v  at up staggered point (i,j+1/2)
          // here we use stale (old) diffusivity if faststep
          if (o==0) {     // If I-offset
            uvbar[o][i][j] = - Df[o][i][j] * slope / thickness + ub[o][i][j];
          } else {        // J-offset
            uvbar[o][i][j] = - Df[o][i][j] * slope / thickness + vb[o][i][j];
          }
         
        } else {  // zero thickness case
          ub[o][i][j] = 0;
          vb[o][i][j] = 0;
          uvbar[o][i][j] = 0;
          if (!faststep) {  // only zero these when would be updated anyway
            Df[o][i][j] = 0;
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
  ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &J[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &J[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[5], &Sigma[1]); CHKERRQ(ierr);

  delete [] delta; delete [] K;
  return 0;
}


//! Compute the rate of frictional heating where the base is sliding by assuming the shear stress is from SIA.
PetscErrorCode IceModel::frictionalHeatingSIAStaggered() {
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub[2], **vb[2], **Rb[2], **mask, **H;

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[8], &Rb[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[9], &Rb[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        // staggered point: o==0 is right, o==1 is up
        const PetscInt      oi = 1-o, oj=o;
        const PetscScalar   thickness = 0.5 * (H[i][j] + H[i+oi][j+oj]);
 
        // basal frictional heating
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          Rb[o][i][j] = 0.0;
        } else { // ignor ice streams; will be overwritten by
                 //   correctBasalFrictionalHeating() if useSSAVelocities==TRUE
          const PetscScalar P = ice.rho * grav * thickness;
          const PetscScalar basal_stress_x = P * h_x[o][i][j];
          const PetscScalar basal_stress_y = P * h_y[o][i][j];
          Rb[o][i][j] = - basal_stress_x * ub[o][i][j] - basal_stress_y * vb[o][i][j];
        }
      }
    }
  }
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &ub[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[5], &ub[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[6], &vb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[7], &vb[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[8], &Rb[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[9], &Rb[1]); CHKERRQ(ierr);
  return 0;
}


//! Put the basal velocity onto the regular grid.
/*!
At the end of velocitySIAStaggered() the basal velocity is available on 
the staggered grid.  This procedure averages it onto the regular grid.  Note that
communication of ghosted values must occur between velocitySIAStaggered() and this 
procedure for the averaging to work.
 */
PetscErrorCode IceModel::basalSIAConditionsToRegular() {
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


//! Put the basal frictional heating and the volume strain-heating onto the regular grid.
/*!
At the end of velocitySIAStaggered() the basal frictional heating and the volume 
strain-heating are available on 
the staggered grid.  This procedure averages them onto the regular grid.  Note that
communication of ghosted values must occur between velocitySIAStaggered() and this 
procedure for the averaging to work.
 */
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


//! Update regular grid horizontal velocities u,v at depth for SIA regions.
/*! 
In the current scheme the procedure velocitySIAStaggered() computes several scalar
quantities at depth (the details of which are too complicated to explain).  That is,
these quantities correspond to three-dimensional arrays.  This procedure takes those 
quantities and computes the three-dimensional arrays for the horizontal components \f$u\f$ and 
\f$v\f$ of the velocity field.  The vertical component \f$w\f$ of the velocity field 
is computed by vertVelocityFromIncompressibility().
 */
PetscErrorCode IceModel::horizontalVelocitySIARegular() {
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

