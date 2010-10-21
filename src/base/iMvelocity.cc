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
  prof->begin(event_sia);

  // do SIA
  if (computeSIAVelocities == PETSC_TRUE) {
    stdout_flags += " SIA ";

    prof->begin(event_sia_2d);
    // uses:
    // * thickness and bed topography or
    // * usurf
    // updates: surface gradient in vWork2d[0,1,2,3], including w=1 ghosts
    ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here
    // surface gradient temporarily stored in vWork2d[0 1 2 3] 
    ierr = verbPrintf(5,grid.com, "{surfaceGradientSIA()}"); CHKERRQ(ierr);

    // uses: 
    // * enthalpy (or temperature) (including w=2 ghosts)
    // * (possibly) age (including w=2 ghosts)
    // * surface gradient (including w=1 ghosts)
    // * thickness (including w=2 ghosts)
    // updates:
    // * SIA velocities on the staggered grid (uvbar), including w=1 ghosts
    // * Sigma on the staggered grid (Sigmastag3[0,1]), including w=1 ghosts
    // * I on the staggered grid (Istag3[0,1])
    // * computes max diffusivity
    ierr = velocitySIAStaggered(); CHKERRQ(ierr);
    ierr = verbPrintf(5,grid.com, "{velocitySIAStaggered()}"); CHKERRQ(ierr);

    // no need to communicate uvbar for boundary conditions for SSA and
    // velocities2DSIAToRegular(): w=1 ghosts were updated locally

    if (mu_sliding == 0.0) { // no need to spend time on nothing
      ierr = vel_basal.set(0.0); CHKERRQ(ierr);
      ierr = vRb.set(0.0); CHKERRQ(ierr);
    } else {
      // compute (and initialize values in) ub,vb and Rb; zero everything where floating
      ierr = basalSlidingHeatingSIA(); CHKERRQ(ierr);
      ierr = verbPrintf(5,grid.com, "{basalSlidingHeatingSIA()}"); CHKERRQ(ierr);
      ierr = vel_basal.beginGhostComm(); CHKERRQ(ierr);
      ierr = vel_basal.endGhostComm(); CHKERRQ(ierr);
      // no need to communicate vRb since not differenced in horizontal
    }

    // now put staggered grid value of vertically-averaged horizontal velocity on regular grid
    // (this makes ubar and vbar be just the result of SIA; note that for SSA
    // we used *saved* ubar,vbar as first guess for iteration at SSA points)
    // also put staggered value of basal velocity onto regular grid

    // uses:
    // * uvbar (including w=1 ghosts)
    // updates:
    // * vel_bar
    ierr = velocities2DSIAToRegular(); CHKERRQ(ierr);
    ierr = verbPrintf(5,grid.com, "{velocities2DSIAToRegular()}"); CHKERRQ(ierr);

    prof->end(event_sia_2d);

    prof->begin(event_sia_3d);
    if (updateVelocityAtDepth) {  
      // I on staggered is used in  horizontalVelocitySIARegular()
      // Sigma on staggered is used in SigmaSIAToRegular()

      // uses:
      // * Sigma on the staggered grid, including w=1 ghosts
      // * thickness
      // updates:
      // * Sigma (local values only; does not need communication)
      ierr = SigmaSIAToRegular(); CHKERRQ(ierr);

      // uses:
      // * surface gradient, including w=1 ghosts
      // * I on the staggered grid, including w=1 ghosts
      // * vel_basal (optional, local values only)
      // updates:
      // * u3, local values only
      // * v3, local values only
      ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
      ierr = verbPrintf(5,grid.com, "{SigmaSIAToRegular(),horizontalVelocitySIARegular()}");
                CHKERRQ(ierr);
    }
    prof->end(event_sia_3d);
  } else { // if computeSIAVelocities == PETSC_FALSE
    ierr = vel_basal.set(0.0); CHKERRQ(ierr);

    // this is a hack; we need to seperate vel_bar, which is a result of the
    // SIA computation, from SSA boundary conditions

    // the following keeps values set as BCs and zeros the rest:
    ierr = vMask.begin_access(); CHKERRQ(ierr);
    ierr = vel_bar.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if (vMask.value(i,j) != MASK_SHEET)
          vel_bar(i,j).u = 0.0;
          vel_bar(i,j).v = 0.0;
      }
    }
    ierr = vel_bar.end_access(); CHKERRQ(ierr);
    ierr = vMask.end_access(); CHKERRQ(ierr);

    ierr = vRb.set(0.0); CHKERRQ(ierr);
    if (updateVelocityAtDepth) {
      ierr = u3.set(0.0); CHKERRQ(ierr);
      ierr = v3.set(0.0); CHKERRQ(ierr);
      ierr = Sigma3.set(0.0); CHKERRQ(ierr);
    }
  }

  prof->end(event_sia);
  PetscLogEventEnd(siaEVENT,0,0,0,0);
  PetscLogEventBegin(ssaEVENT,0,0,0,0);
  prof->begin(event_ssa);

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

    // uses:
    // * mask
    // * vel_ssa
    // * vel_basal
    // * vel_bar
    // * u3
    // * v3
    // * uvbar, including w=1 ghosts
    // updates:
    // * u3
    // * v3
    // * vel_bar
    // * vel_basal
    // (In all these cases only local values are updated.)
    ierr = broadcastSSAVelocity(updateVelocityAtDepth); CHKERRQ(ierr);

    // now communicate modified velocity fields
    ierr = vel_basal.beginGhostComm(); CHKERRQ(ierr);
    ierr = vel_basal.endGhostComm(); CHKERRQ(ierr);

    // note correctSigma() differences ubasal,vbasal in horizontal, so communication
    //   above is important
    ierr = correctSigma(); CHKERRQ(ierr);
    ierr = correctBasalFrictionalHeating(); CHKERRQ(ierr);
  } // end of if (use_ssa_velocity) { ...

  prof->end(event_ssa);
PetscLogEventEnd(ssaEVENT,0,0,0,0);
PetscLogEventBegin(velmiscEVENT,0,0,0,0);

  prof->begin(event_vel_com);
  // in latter case u,v are modified by broadcastSSAVelocity():
  if (updateVelocityAtDepth) {  
    ierr = u3.beginGhostComm(); CHKERRQ(ierr);
    ierr = v3.beginGhostComm(); CHKERRQ(ierr);
    ierr = u3.endGhostComm(); CHKERRQ(ierr);
    ierr = v3.endGhostComm(); CHKERRQ(ierr);
  }
  prof->end(event_vel_com);

  prof->begin(event_vel_inc);
  // finally update w
  if (updateVelocityAtDepth) {
    ierr = vertVelocityFromIncompressibility(); CHKERRQ(ierr);
    // no communication needed for w, which is only differenced in the column
  }
  prof->end(event_vel_inc);
  
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
location of the base of the ice column</i>.  That is, the vertical velocity
computed here is identified as \f$\tilde w(x,y,s,t)\f$ in the page
\ref vertchange.

Thus \f$w<0\f$ here means that that
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
  ierr = vbmr.begin_access(); CHKERRQ(ierr);

  PetscScalar *w, *u_im1, *u_ip1, *v_jm1, *v_jp1;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = w3.getInternalColumn(i,j,&w); CHKERRQ(ierr);

      ierr = u3.getInternalColumn(i-1,j,&u_im1); CHKERRQ(ierr);
      ierr = u3.getInternalColumn(i+1,j,&u_ip1); CHKERRQ(ierr);

      ierr = v3.getInternalColumn(i,j-1,&v_jm1); CHKERRQ(ierr);
      ierr = v3.getInternalColumn(i,j+1,&v_jp1); CHKERRQ(ierr);

      if (include_bmr_in_continuity) {
        w[0] = - vbmr(i,j);
      } else {
        w[0] = 0.0;
      }

      PetscScalar OLDintegrand
             = (u_ip1[0] - u_im1[0]) / (2.0*dx) + (v_jp1[0] - v_jm1[0]) / (2.0*dy);
      for (PetscInt k = 1; k < grid.Mz; ++k) {
        const PetscScalar NEWintegrand
             = (u_ip1[k] - u_im1[k]) / (2.0*dx) + (v_jp1[k] - v_jm1[k]) / (2.0*dy);
        const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
        w[k] = w[k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;
      }
    }
  }

  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  return 0;
}


