// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"
#include <petscvec.h>


//! Compute the maximum diffusivity associated to the SIA deformational velocity.
/*! 
The time-stepping scheme for mass continuity is explicit.

For the non-sliding, deformational part of the vertically-integrated 
horizontal mass flux \f$\mathbf{q}\f$, the partial differential equation 
is diffusive.  Thus there is a stability criterion \ref MortonMayers 
which depends on the diffusivity coefficient.  Of course,
because the PDE is nonlinear, this diffusivity changes at every time step.  This 
procedure computes the maximum of the diffusivity on the grid.

See determineTimeStep() and massContExplicitStep().
 */
PetscErrorCode IceModel::computeMaxDiffusivity(bool updateDiffusViewer) {
  // assumes vuvbar holds correct deformational values of velocities

  PetscErrorCode ierr;

  const PetscScalar DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE = 1.0e-4;
  PetscScalar **h, **H, **uvbar[2], **D;
  PetscScalar Dmax = 0.0;

  ierr =        vh.get_array(h); CHKERRQ(ierr);
  ierr =        vH.get_array(H); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(D); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        if (computeSIAVelocities == PETSC_TRUE) {
          // note: when basal sliding is proportional to surface slope, as
          // it usually will be when sliding occurs in a MASK_SHEET area, then
          //    D = H Ubar / alpha
          // is the correct formula; note division by zero is avoided by
          // addition to alpha
          const PetscScalar h_x=(h[i+1][j]-h[i-1][j])/(2.0*grid.dx),
                            h_y=(h[i][j+1]-h[i][j-1])/(2.0*grid.dy),
                            alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
          const PetscScalar udef = 0.5 * (uvbar[0][i][j] + uvbar[0][i-1][j]),
                            vdef = 0.5 * (uvbar[1][i][j] + uvbar[1][i][j-1]),
                            Ubarmag = sqrt(PetscSqr(udef) + PetscSqr(vdef));
          const PetscScalar d =
               H[i][j] * Ubarmag/(alpha + DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE);
          if (d > Dmax) Dmax = d;
          D[i][j] = d;
        } else {
          D[i][j] = 0.0; // no diffusivity if no SIA
        }
      } else {
        D[i][j] = 0.0; // no diffusivity if no ice
      }
    }
  }
  ierr =        vh.end_access(); CHKERRQ(ierr);
  ierr =        vH.end_access(); CHKERRQ(ierr);  
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  if (updateDiffusViewer) { // view diffusivity (m^2/s)
    ierr = update2DViewer('D',vWork2d[0],1.0); CHKERRQ(ierr);
  }

  ierr = PetscGlobalMax(&Dmax, &gDmax, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Compute the maximum velocities for time-stepping and reporting to user.
/*!
Computes the maximum magnitude of the components \f$u,v,w\f$ of the 3D velocity.
Then sets \c CFLmaxdt, the maximum time step allowed under the 
Courant-Friedrichs-Lewy (CFL) condition on the 
horizontal advection scheme for age and for temperature.

Under BOMBPROOF there is no CFL condition for the vertical advection.
The maximum vertical velocity is computed but it does not affect
\c CFLmaxdt.
 */
PetscErrorCode IceModel::computeMax3DVelocities() {
  PetscErrorCode ierr;
  PetscScalar **H, *u, *v, *w;
  PetscScalar locCFLmaxdt = maxdt;

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);

  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt      ks = grid.kBelowHeight(H[i][j]);
      ierr = u3.getInternalColumn(i,j,&u); CHKERRQ(ierr);
      ierr = v3.getInternalColumn(i,j,&v); CHKERRQ(ierr);
      ierr = w3.getInternalColumn(i,j,&w); CHKERRQ(ierr);
      for (PetscInt k=0; k<ks; ++k) {
        const PetscScalar absu = PetscAbs(u[k]),
                          absv = PetscAbs(v[k]);
        maxu = PetscMax(maxu,absu);
        maxv = PetscMax(maxv,absv);
        // make sure the denominator below is positive:
        PetscScalar tempdenom = (0.001/secpera) / (grid.dx + grid.dy);  
        tempdenom += PetscAbs(absu/grid.dx) + PetscAbs(absv/grid.dy);
        locCFLmaxdt = PetscMin(locCFLmaxdt,1.0 / tempdenom); 
        maxw = PetscMax(maxw, PetscAbs(w[k]));        
      }
    }
  }

  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxu, &gmaxu, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxv, &gmaxv, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxw, &gmaxw, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMin(&locCFLmaxdt, &CFLmaxdt, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Compute the CFL constant associated to first-order upwinding for the sliding contribution to mass continuity.
/*!
This procedure computes the maximum horizontal speed in the SSA areas.  In
particular it computes CFL constant for the upwinding, in massContExplicitStep(),
which applies to the basal component of mass flux.

That is, because the map-plane mass continuity is advective in the
sliding case we have a CFL condition.
 */
PetscErrorCode IceModel::computeMax2DSlidingSpeed() {
  PetscErrorCode ierr;
  PetscScalar **ub, **vb, **mask;
  PetscScalar locCFLmaxdt2D = maxdt;
  
  ierr = vub.get_array(ub); CHKERRQ(ierr);
  ierr = vvb.get_array(vb); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // the following conditionals, both -ocean_kill and -float_kill, are also applied in 
      //   IceModel::massContExplicitStep() when zeroing thickness
      const bool ignorableOcean =
          (   ( (doOceanKill == PETSC_TRUE) && (PismIntMask(mask[i][j]) == MASK_FLOATING_OCEAN0) )
           || ( (floatingIceKilled == PETSC_TRUE) && (PismModMask(mask[i][j]) == MASK_FLOATING)  )  );
      if (!ignorableOcean) {
        PetscScalar denom = PetscAbs(ub[i][j])/grid.dx + PetscAbs(vb[i][j])/grid.dy;
        denom += (0.01/secpera)/(grid.dx + grid.dy);  // make sure it's pos.
        locCFLmaxdt2D = PetscMin(locCFLmaxdt2D,1.0/denom);
      }
    }
  }
  ierr = vub.end_access(); CHKERRQ(ierr);
  ierr = vvb.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalMin(&locCFLmaxdt2D, &CFLmaxdt2D, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Compute the maximum time step allowed by the diffusive SIA.
/*!
Note computeMaxDiffusivity() must be called before this to set \c gDmax.  Note
adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio * dx^2/maxD
(if dx=dy)

Reference: \ref MortonMayers pp 62--63.
 */
PetscErrorCode IceModel::adaptTimeStepDiffusivity() {
  const PetscScalar DEFAULT_ADDED_TO_GDMAX_ADAPT = 1.0e-2;
  const PetscScalar  
          gridfactor = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dt_from_diffus = adaptTimeStepRatio
                     * 2 / ((gDmax + DEFAULT_ADDED_TO_GDMAX_ADAPT) * gridfactor);
  if ((doSkip == PETSC_TRUE) && (skipCountDown == 0)) {
    const PetscScalar  conservativeFactor = 0.95;
    // typically "dt" in next line is from CFL for advection in temperature equation,
    //   but in fact it might be from other restrictions, e.g. CFL for mass continuity
    //   in basal sliding case, or max_dt
    skipCountDown = (PetscInt) floor(conservativeFactor * (dt / dt_from_diffus));
    skipCountDown = ( skipCountDown >  skipMax) ?  skipMax :  skipCountDown;
  } // if  skipCountDown > 0 then it will get decremented at the mass balance step
  if (dt_from_diffus < dt) {
    dt = dt_from_diffus;
    adaptReasonFlag = 'd';
  }
  return 0;
}


//! Use various stability criteria to determine the time step for an evolution run.
/*! 
The main loop in run() approximates many physical processes.  Several of these approximations,
including the mass continuity and temperature equations in particular, involve stability
criteria.  This procedure builds the length of the next time step by using these criteria and 
by incorporating choices made by options (e.g. <c>-max_dt</c>) and by derived classes.
 */
PetscErrorCode IceModel::determineTimeStep(const bool doTemperatureCFL) {
  PetscErrorCode ierr;

  if ( (runtimeViewers[cIndex('D')] != PETSC_NULL) 
       || ( (doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE) ) ) {
    ierr = computeMaxDiffusivity(true); CHKERRQ(ierr);
  }
  const PetscScalar timeToEnd = (endYear-grid.year) * secpera;
  if (dt_force > 0.0) {
    dt = dt_force; // override usual dt mechanism
    adaptReasonFlag = 'f';
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
  } else {
    dt = maxdt;
    adaptReasonFlag = 'm';
    if (doAdaptTimeStep == PETSC_TRUE) {
      if ((doTemp == PETSC_TRUE) && (doTemperatureCFL == PETSC_TRUE)) {
        // CFLmaxdt is set by computeMax3DVelocities() in call to velocity() iMvelocity.cc
        dt_from_cfl = CFLmaxdt;
        if (dt_from_cfl < dt) {
          dt = dt_from_cfl;
          adaptReasonFlag = 'c';
        }
      }
      if ((doMassConserve == PETSC_TRUE) && (useSSAVelocity == PETSC_TRUE)) {
        // CFLmaxdt2D is set by broadcastSSAVelocity()
        if (CFLmaxdt2D < dt) {
          dt = CFLmaxdt2D;
          adaptReasonFlag = 'u';
        }
      }
      if ((doMassConserve == PETSC_TRUE) && (computeSIAVelocities == PETSC_TRUE)) {
        // note: if doSkip then skipCountDown = floor(dt_from_cfl/dt_from_diffus)
        ierr = adaptTimeStepDiffusivity(); CHKERRQ(ierr); // might set adaptReasonFlag = 'd'
      }
    }
    if ((maxdt_temporary > 0.0) && (maxdt_temporary < dt)) {
      dt = maxdt_temporary;
      adaptReasonFlag = 't';
    }
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
    if ((adaptReasonFlag == 'm') || (adaptReasonFlag == 't') || (adaptReasonFlag == 'e')) {
      if (skipCountDown > 1) skipCountDown = 1; 
    }
  }    
  return 0;
}

