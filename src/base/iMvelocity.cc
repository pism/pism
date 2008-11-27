// Copyright (C) 2004-2008 Jed Brown, Ed Bueler and Constantine Khroulev
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


//! Manage the computation of velocity and do the necessary parallel communication.
/*!
This procedure calls the important routines velocitySIAStaggered(), velocitySSA(),
and vertVelocityFromIncompressibility() according to various flags. 
 */
PetscErrorCode IceModel::velocity(bool updateVelocityAtDepth) {
  PetscErrorCode ierr;
  static PetscTruth firstTime = PETSC_TRUE;

#if (PISM_LOG_EVENTS)
PetscLogEventBegin(siaEVENT,0,0,0,0);
#endif

  // do SIA
  if (computeSIAVelocities == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com, " SIA "); CHKERRQ(ierr);
    ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here for eta
    // surface gradient temporarily stored in vWork2d[0 1 2 3] 
    ierr = verbPrintf(5,grid.com, "{surfaceGradientSIA()}"); CHKERRQ(ierr);

    // communicate h_x[o], h_y[o] on staggered for basalSIA() and horizontalVelocitySIARegular()
    for (PetscInt k=0; k<4; ++k) { 
      ierr = vWork2d[k].beginGhostComm(); CHKERRQ(ierr);
      ierr = vWork2d[k].endGhostComm(); CHKERRQ(ierr);
    }

    ierr = velocitySIAStaggered(); CHKERRQ(ierr);
    ierr = verbPrintf(5,grid.com, "{velocitySIAStaggered()}"); CHKERRQ(ierr);

    // communicate vuvbar[01] for boundary conditions for SSA and vertAveragedVelocityToRegular()
    // and velocities2DSIAToRegular()
    ierr = vuvbar[0].beginGhostComm();
    ierr = vuvbar[0].endGhostComm();
    ierr = vuvbar[1].beginGhostComm();
    ierr = vuvbar[1].endGhostComm();
    ierr = verbPrintf(5,grid.com, "{comm after velocitySIAStaggered()}"); CHKERRQ(ierr);

    if (muSliding == 0.0) { // no need to spend time on nothing
      ierr = vub.set(0.0); CHKERRQ(ierr);
      ierr = vvb.set(0.0); CHKERRQ(ierr);
      ierr = vRb.set(0.0); CHKERRQ(ierr);
    } else {
      // compute (and initialize values in) ub,vb and Rb; zero everything where floating
      ierr = basalSlidingHeatingSIA(); CHKERRQ(ierr);
      ierr = verbPrintf(5,grid.com, "{basalSlidingHeatingSIA()}"); CHKERRQ(ierr);
      ierr = vub.beginGhostComm(); CHKERRQ(ierr);
      ierr = vub.endGhostComm(); CHKERRQ(ierr);
      ierr = vvb.beginGhostComm(); CHKERRQ(ierr);
      ierr = vvb.endGhostComm(); CHKERRQ(ierr);
      // no need to communicate vRb since not differenced in horizontal
    }

    // now put staggered grid value of vertically-averaged horizontal velocity on regular grid
    // (this makes ubar and vbar be just the result of SIA; note that for SSA
    // we used *saved* ubar,vbar as first guess for iteration at SSA points)
    // also put staggered value of basal velocity onto regular grid
    ierr = velocities2DSIAToRegular(); CHKERRQ(ierr);
    ierr = verbPrintf(5,grid.com, "{velocities2DSIAToRegular()}"); CHKERRQ(ierr);
    ierr = vubar.beginGhostComm(); CHKERRQ(ierr);
    ierr = vubar.endGhostComm(); CHKERRQ(ierr);
    ierr = vvbar.beginGhostComm(); CHKERRQ(ierr);
    ierr = vvbar.endGhostComm(); CHKERRQ(ierr);
  
    if (updateVelocityAtDepth) {  
      // communicate I on staggered for horizontalVelocitySIARegular()
      //   and also communicate Sigma on staggered for SigmaSIAToRegular()
      ierr = Sigmastag3[0].beginGhostComm(); CHKERRQ(ierr);
      ierr = Sigmastag3[1].beginGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
      ierr = Sigmastag3[0].endGhostComm(); CHKERRQ(ierr);
      ierr = Sigmastag3[1].endGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);

      ierr = SigmaSIAToRegular(); CHKERRQ(ierr);
      ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
      ierr = verbPrintf(5,grid.com, "{SigmaSIAToRegular(),horizontalVelocitySIARegular()}");
                CHKERRQ(ierr);
    }
  } else { // if computeSIAVelocities == PETSC_FALSE
    // do NOT zero out vuvbar[0],vuvbar[1]; they are used to communicate boundary
    // conditions to SSA calculation
    ierr = vubar.set(0.0); CHKERRQ(ierr);
    ierr = vvbar.set(0.0); CHKERRQ(ierr);
    ierr = vub.set(0.0); CHKERRQ(ierr);
    ierr = vvb.set(0.0); CHKERRQ(ierr);
    ierr = vRb.set(0.0); CHKERRQ(ierr);
    if (updateVelocityAtDepth) {
      ierr = u3.set(0.0); CHKERRQ(ierr);
      ierr = v3.set(0.0); CHKERRQ(ierr);
      ierr = Sigma3.set(0.0); CHKERRQ(ierr);
    }
    ierr = verbPrintf(2,grid.com, "     "); CHKERRQ(ierr);
  }

#if (PISM_LOG_EVENTS)
PetscLogEventEnd(siaEVENT,0,0,0,0);
PetscLogEventBegin(ssaEVENT,0,0,0,0);
#endif
  
  // do SSA
  if ((firstTime == PETSC_TRUE) && (useSSAVelocity == PETSC_TRUE)) {
    ierr = initSSA(); CHKERRQ(ierr);
  }
  if (useSSAVelocity == PETSC_TRUE) { // communication happens within SSA
    if ((firstTime == PETSC_TRUE) || (updateVelocityAtDepth)) {
      ierr = verbPrintf(2,grid.com, "SSA"); CHKERRQ(ierr);
      PetscInt numSSAiter;
      ierr = velocitySSA(&numSSAiter); CHKERRQ(ierr); // comm here ...
      ierr = verbPrintf(2,grid.com," "); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com,"       "); CHKERRQ(ierr);
    }
    // even if velocitySSA() did not run, we still need to use stored SSA velocities 
    // to get 3D velocity field, basal velocities, basal frictional heating, 
    // and strain dissipation heating
    ierr = broadcastSSAVelocity(updateVelocityAtDepth); CHKERRQ(ierr);

    // now communicate modified velocity fields
    ierr = vubar.beginGhostComm(); CHKERRQ(ierr);
    ierr = vubar.endGhostComm(); CHKERRQ(ierr);
    ierr = vvbar.beginGhostComm(); CHKERRQ(ierr);
    ierr = vvbar.endGhostComm(); CHKERRQ(ierr);

    ierr = vub.beginGhostComm(); CHKERRQ(ierr);
    ierr = vub.endGhostComm(); CHKERRQ(ierr);
    ierr = vvb.beginGhostComm(); CHKERRQ(ierr);
    ierr = vvb.endGhostComm(); CHKERRQ(ierr);

    // note correctSigma() differences ub,vb in horizontal, so communication
    //   above is important
    ierr = correctSigma(); CHKERRQ(ierr);
    ierr = correctBasalFrictionalHeating(); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, "       "); CHKERRQ(ierr);
  }

#if (PISM_LOG_EVENTS)
PetscLogEventEnd(ssaEVENT,0,0,0,0);
PetscLogEventBegin(velmiscEVENT,0,0,0,0);
#endif

  // in latter case u,v are modified by broadcastSSAVelocity():
  if (updateVelocityAtDepth) {  
    ierr = u3.beginGhostComm(); CHKERRQ(ierr);
    ierr = v3.beginGhostComm(); CHKERRQ(ierr);
    ierr = u3.endGhostComm(); CHKERRQ(ierr);
    ierr = v3.endGhostComm(); CHKERRQ(ierr);
  }

  // finally update w
  if (updateVelocityAtDepth) {
    ierr = vertVelocityFromIncompressibility(); CHKERRQ(ierr);
    // no communication needed for w, which is only differenced in the column
  }
  
  // communication here for global max; sets CFLmaxdt2D
  ierr = computeMax2DSlidingSpeed(); CHKERRQ(ierr);   

  if (updateVelocityAtDepth) {
    // communication here for global max; sets CFLmaxdt
    ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  }

#if (PISM_LOG_EVENTS)
PetscLogEventEnd(velmiscEVENT,0,0,0,0);
#endif
  
  firstTime = PETSC_FALSE;
  return 0;
}


//! Compute vertical velocity using basal conditions and incompressibility of the ice.
/*! 
The original statement of incompressibility is
   \f[ \nabla\cdot\mathbf{U} + \frac{\partial w}{\partial z} = 0. \f]
This is immediately equivalent to the integral
   \f[ w(x,y,z,t) = - \int_{b(x,y,t)}^{z} \nabla\cdot\mathbf{U}\,d\zeta
                           + w_b(x,y,t). \f]

The basal kinematic equation is
   \f[ w_b = \frac{\partial b}{\partial t} + \mathbf{U}_b \cdot \nabla b - S \f]
where \f$S\f$ is the basal melt rate.  (The inclusion of the basal melt rate is optional.)  
This equation determines the vertical velocity of the ice at the base (\f$w_b\f$), 
which is needed in the incompressibility integral.

@cond VERTCHANGE
Note we do a change of vertical coordinate \f$\tilde z = z - b(x,y,t)\f$.  Thus all 
derivatives, with respect to any variable \f$x,y,z,t\f$, change accordingly 
(see elsewhere).  The revised form of the incompressibility integral is
   \f[ w(x,y,\tilde z,t) = - \int_0^{\tilde z} \left(
          \nabla\cdot\mathbf{U} 
          - \nabla b \cdot \left(\frac{\partial u}{\partial \tilde z},
                                 \frac{\partial v}{\partial \tilde z}\right)
                                                     \right)\,d\zeta
                           + w_b(x,y,t). \f]
@endcond

The vertical integral is computed by the trapezoid rule.  There is no assumption about equal
spacing.
 */
PetscErrorCode IceModel::vertVelocityFromIncompressibility() {
  PetscErrorCode  ierr;
  const PetscScalar   dx = grid.dx, 
                      dy = grid.dy;
  const PetscInt      Mz = grid.Mz;
  PetscScalar **ub, **vb, **basalMeltRate, **mask, **b, **dbdt;
  PetscScalar *u, *v, *w;

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);

  ierr = vub.get_array(ub); CHKERRQ(ierr);
  ierr = vvb.get_array(vb); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vbed.get_array(b); CHKERRQ(ierr);
  ierr = vuplift.get_array(dbdt); CHKERRQ(ierr);
  ierr = vbasalMeltRate.get_array(basalMeltRate); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = u3.getInternalColumn(i,j,&u); CHKERRQ(ierr);
      ierr = v3.getInternalColumn(i,j,&v); CHKERRQ(ierr);
      ierr = w3.getInternalColumn(i,j,&w); CHKERRQ(ierr);

      PetscScalar dbdx = 0.0, dbdy = 0.0;
      if (modMask(mask[i][j]) != MASK_FLOATING) {
        // bed gradient needed if grounded; there is no reason to assume that the
        //   bed elevation has been periodized by the user
        if (i == 0) {
          dbdx = (b[i+1][j] - b[i][j]) / (dx);
        } else if (i == grid.Mx-1) {
          dbdx = (b[i][j] - b[i-1][j]) / (dx);
        } else {
          dbdx = (b[i+1][j] - b[i-1][j]) / (2.0*dx);
        }
        if (j == 0) {
          dbdy = (b[i][j+1] - b[i][j]) / (dy);
        } else if (j == grid.My-1) {
          dbdy = (b[i][j] - b[i][j-1]) / (dy);
        } else {
          dbdy = (b[i][j+1] - b[i][j-1]) / (2.0*dy);
        }
      }

      // basal w from basal kinematical equation
//      w[0] = dbdt[i][j] + ub[i][j] * dbdx + vb[i][j] * dbdy; // DEBUG?: remove dbdt
      w[0] = ub[i][j] * dbdx + vb[i][j] * dbdy;

      if (includeBMRinContinuity == PETSC_TRUE) {
        w[0] -= basalMeltRate[i][j];
      }

      // compute w above base by trapezoid rule 
      planeStar uss, vss;
      ierr = u3.getPlaneStarZ(i,j,0.0,&uss);
      ierr = v3.getPlaneStarZ(i,j,0.0,&vss);
      PetscScalar OLDintegrand = (uss.ip1 - uss.im1) / (2.0*dx) + (vss.jp1 - vss.jm1) / (2.0*dy);
      // at bottom, difference up:
      PetscScalar dz = grid.zlevels[1] - grid.zlevels[0];
      OLDintegrand -= dbdx * (u[1] - u[0]) / dz + dbdy * (v[1] - v[0]) / dz;

      for (PetscInt k = 1; k < Mz; ++k) {
        ierr = u3.getPlaneStarZ(i,j,grid.zlevels[k],&uss);
        ierr = v3.getPlaneStarZ(i,j,grid.zlevels[k],&vss);
        PetscScalar NEWintegrand = (uss.ip1 - uss.im1) / (2.0*dx) + (vss.jp1 - vss.jm1) / (2.0*dy);
        dz = grid.zlevels[k] - grid.zlevels[k-1];
        if (k == Mz-1) { // at top, difference down:
          NEWintegrand -= dbdx * (u[k] - u[k-1]) / dz + dbdy * (v[k] - v[k-1]) / dz;
        } else { // usual case; central difference
          const PetscScalar twodz = grid.zlevels[k+1] - grid.zlevels[k-1];
          NEWintegrand -= dbdx * (u[k+1] - u[k-1]) / twodz + dbdy * (v[k+1] - v[k-1]) / twodz;
        }
        w[k] = w[k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;
      }
      // no need to call w3.setInternalColumn; already set
    }
  }

  ierr =     vub.end_access(); CHKERRQ(ierr);
  ierr =     vvb.end_access(); CHKERRQ(ierr);
  ierr =   vMask.end_access(); CHKERRQ(ierr);
  ierr =    vbed.end_access(); CHKERRQ(ierr);
  ierr = vuplift.end_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);
  
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  return 0;
}


