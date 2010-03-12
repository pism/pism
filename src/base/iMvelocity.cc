// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

  double mu_sliding = config.get("mu_sliding");

  bool use_ssa_velocity = config.get_flag("use_ssa_velocity");

PetscLogEventBegin(siaEVENT,0,0,0,0);

  // do SIA
  if (computeSIAVelocities == PETSC_TRUE) {
    stdout_flags += " SIA ";
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

    if (mu_sliding == 0.0) { // no need to spend time on nothing
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
  }

PetscLogEventEnd(siaEVENT,0,0,0,0);
PetscLogEventBegin(ssaEVENT,0,0,0,0);
  
  // do SSA
  if ((firstTime == PETSC_TRUE) && use_ssa_velocity) {
    ierr = initSSA(); CHKERRQ(ierr);
  }
  if (use_ssa_velocity) { // communication happens within SSA
    if ((firstTime == PETSC_TRUE) || (updateVelocityAtDepth)) {
      PetscInt numSSAiter;
      ierr = velocitySSA(&numSSAiter); CHKERRQ(ierr); // comm here ...
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
  }

PetscLogEventEnd(ssaEVENT,0,0,0,0);
PetscLogEventBegin(velmiscEVENT,0,0,0,0);

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

PetscLogEventEnd(velmiscEVENT,0,0,0,0);
  
  firstTime = PETSC_FALSE;
  return 0;
}


//! Compute vertical velocity using incompressibility of the ice.
/*!
The vertical velocity \f$w(x,y,z,t)\f$ is the velocity <i>relative to the
location of the base of the ice column</i>.  Thus \f$w<0\f$ means that that
that part of the ice is getting closer to the base of the ice, and so on.
The slope of the bed (i.e. relative to the geoid) and/or the motion of the
bed (i.e. from bed deformation) do not affect the vertical velocity.

In fact the following statement is exactly true if the basal melt rate is zero:
the vertical velocity at a point in the ice is positive (negative) if and only
if the average horizontal divergence of the horizontal velocity, in the portion
of the ice column below that point, is negative (positive).
In particular, because \f$z=0\f$ is the location of the base of the ice
always, the only way to have \f$w(x,y,0,t) \ne 0\f$ is to have a basal melt
rate.

Incompressibility itself says
   \f[ \nabla\cdot\mathbf{U} + \frac{\partial w}{\partial z} = 0. \f]
This is immediately equivalent to the integral
   \f[ w(x,y,z,t) = - \int_{b(x,y,t)}^{z} \nabla\cdot\mathbf{U}\,d\zeta
                           + w_b(x,y,t). \f]
Here the value \f$w_b(x,y,t)\f$ is either zero or the negative of the basal melt rate
according to the value of the flag \c include_bmr_in_continuity.

The vertical integral is computed by the trapezoid rule.
 */
PetscErrorCode IceModel::vertVelocityFromIncompressibility() {
  PetscErrorCode  ierr;
  const PetscScalar dx = grid.dx, dy = grid.dy;
  bool include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity");

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.begin_access(); CHKERRQ(ierr);

  PetscScalar *w;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = w3.getInternalColumn(i,j,&w); CHKERRQ(ierr);
      if (include_bmr_in_continuity) {
        w[0] = - vbasalMeltRate(i,j);
      } else {
        w[0] = 0.0;
      }

      planeStar uss, vss;
      ierr = u3.getPlaneStarZ(i,j,0.0,&uss);
      ierr = v3.getPlaneStarZ(i,j,0.0,&vss);
      PetscScalar OLDintegrand
             = (uss.ip1 - uss.im1) / (2.0*dx) + (vss.jp1 - vss.jm1) / (2.0*dy);
      for (PetscInt k = 1; k < grid.Mz; ++k) {
        ierr = u3.getPlaneStarZ(i,j,grid.zlevels[k],&uss);
        ierr = v3.getPlaneStarZ(i,j,grid.zlevels[k],&vss);
        const PetscScalar NEWintegrand
             = (uss.ip1 - uss.im1) / (2.0*dx) + (vss.jp1 - vss.jm1) / (2.0*dy);
        const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
        w[k] = w[k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;
      }
    }
  }

  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  return 0;
}


