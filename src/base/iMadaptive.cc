// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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
is diffusive.  Thus there is a stability criterion \lo\cite{MortonMayers}\elo 
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

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);
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
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);

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

/*
  // compute dzEQ which will be used inside temperatureStep() and ageStep()
  PetscInt    Mz_for_dzEQ, dummyM;
  ierr = getMzMbzForTempAge(Mz_for_dzEQ,dummyM); CHKERRQ(ierr);
  const PetscScalar dzEQ = grid.Lz / ((PetscScalar) (Mz_for_dzEQ - 1));
*/

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);

  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt      ks = grid.kBelowHeight(H[i][j]);
/*
      const bool isMarginal = checkThinNeigh(
             H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
*/
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
/*
        if (!isMarginal) {
          const PetscScalar absw = PetscAbs(w[k]);
          maxw = PetscMax(maxw,absw);
          tempdenom += PetscAbs(absw / dzEQ);
        }
*/
        locCFLmaxdt = PetscMin(locCFLmaxdt,1.0 / tempdenom); 
        maxw = PetscMax(maxw, PetscAbs(w[k]));        
      }
    }
  }

  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = w3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

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
  PetscScalar **ub, **vb;
  PetscScalar locCFLmaxdt2D = maxdt;
  
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar denom = PetscAbs(ub[i][j])/grid.dx + PetscAbs(vb[i][j])/grid.dy;
      denom += (0.01/secpera)/(grid.dx + grid.dy);  // make sure it's pos.
      locCFLmaxdt2D = PetscMin(locCFLmaxdt2D,1.0/denom);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);

  ierr = PetscGlobalMin(&locCFLmaxdt2D, &CFLmaxdt2D, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Compute the maximum time step allowed by the diffusive SIA.
/*!
Note computeMaxDiffusivity() must be called before this to set \c gDmax.  Note
adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio * dx^2/maxD
(if dx=dy)

Reference: \lo\cite{MortonMayers}\elo pp 62--63.
 */
PetscErrorCode IceModel::adaptTimeStepDiffusivity() {
  const PetscScalar DEFAULT_ADDED_TO_GDMAX_ADAPT = 1.0e-2;
  const PetscScalar  
          gridfactor = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dt_from_diffus = adaptTimeStepRatio
                     * 2 / ((gDmax + DEFAULT_ADDED_TO_GDMAX_ADAPT) * gridfactor);
  if ((doTempSkip == PETSC_TRUE) && (tempskipCountDown == 0)) {
//    const PetscScalar  conservativeFactor = 0.8;
    const PetscScalar  conservativeFactor = 0.95;
    // typically "dt" in next line is from CFL, but might be from other, e.g. maxdt
    tempskipCountDown = (PetscInt) floor(conservativeFactor * (dt / dt_from_diffus));
    tempskipCountDown = (tempskipCountDown > tempskipMax) ? tempskipMax : tempskipCountDown;
  } // if tempskipCountDown > 0 then it will get decremented at the mass balance step
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
    if ((doAdaptTimeStep == PETSC_TRUE) && (doTemp == PETSC_TRUE)
        && doTemperatureCFL) {
      // CFLmaxdt is set by computeMax3DVelocities() in call to velocity() iMvelocity.cc
      dt_from_cfl = CFLmaxdt;
      if (dt_from_cfl < dt) {
        dt = dt_from_cfl;
        adaptReasonFlag = 'c';
      }
    } 
    if ((doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE)
        && (useSSAVelocity)) {
      // CFLmaxdt2D is set by broadcastSSAVelocity()
      if (CFLmaxdt2D < dt) {
        dt = CFLmaxdt2D;
        adaptReasonFlag = 'u';
      }
    }
    if ((doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE)) {
      // note: if doTempSkip then tempskipCountDown = floor(dt_from_cfl/dt_from_diffus)
      ierr = adaptTimeStepDiffusivity(); CHKERRQ(ierr); // might set adaptReasonFlag = 'd'
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
      if (tempskipCountDown > 1) tempskipCountDown = 1; 
    }
  }    
  return 0;
}

