// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "bedrockThermalUnit.hh"
#include "PISMTime.hh"

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
  PetscScalar *u, *v, *w;
  PetscScalar locCFLmaxdt = config.get("maximum_time_step_years", "years", "seconds");

  IceModelVec3 *u3, *v3, *w3;

  MaskQuery mask(vMask);

  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      if (mask.icy(i, j)) {
        const PetscInt ks = grid.kBelowHeight(vH(i, j));
        ierr = u3->getInternalColumn(i, j, &u); CHKERRQ(ierr);
        ierr = v3->getInternalColumn(i, j, &v); CHKERRQ(ierr);
        ierr = w3->getInternalColumn(i, j, &w); CHKERRQ(ierr);
        for (PetscInt k = 0; k <= ks; ++k) {
          const PetscScalar absu = PetscAbs(u[k]),
            absv = PetscAbs(v[k]);
          maxu = PetscMax(maxu, absu);
          maxv = PetscMax(maxv, absv);
          // make sure the denominator below is positive:
          PetscScalar tempdenom = (0.001 / secpera) / (grid.dx + grid.dy);
          tempdenom += PetscAbs(absu / grid.dx) + PetscAbs(absv / grid.dy);
          locCFLmaxdt = PetscMin(locCFLmaxdt, 1.0 / tempdenom);
          maxw = PetscMax(maxw, PetscAbs(w[k]));
        }
      }
    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&maxu, &gmaxu, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxv, &gmaxv, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxw, &gmaxw, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMin(&locCFLmaxdt, &CFLmaxdt, grid.com); CHKERRQ(ierr);
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
  PISMVector2 **vel;
  PetscScalar locCFLmaxdt2D = config.get("maximum_time_step_years",
                                         "years", "seconds");

  MaskQuery mask(vMask);

  IceModelVec2V *vel_advective;
  ierr = stress_balance->get_2D_advective_velocity(vel_advective); CHKERRQ(ierr);

  ierr = vel_advective->get_array(vel); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.icy(i, j)) {
        PetscScalar denom = PetscAbs(vel[i][j].u)/grid.dx + PetscAbs(vel[i][j].v)/grid.dy;
        denom += (0.01/secpera)/(grid.dx + grid.dy);  // make sure it's pos.
        locCFLmaxdt2D = PetscMin(locCFLmaxdt2D,1.0/denom);
      }
    }
  }
  ierr = vel_advective->end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMin(&locCFLmaxdt2D, &CFLmaxdt2D, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Compute the maximum time step allowed by the diffusive SIA.
/*! Note adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio *
dx^2/maxD (if dx=dy)

Reference: [\ref MortonMayers] pp 62--63.
 */
PetscErrorCode IceModel::adaptTimeStepDiffusivity() {
  PetscErrorCode ierr;

  bool do_skip = config.get_flag("do_skip");

  const PetscScalar adaptTimeStepRatio = config.get("adaptive_timestepping_ratio");

  const PetscScalar DEFAULT_ADDED_TO_GDMAX_ADAPT = 1.0e-2;

  ierr = stress_balance->get_max_diffusivity(gDmax); CHKERRQ(ierr);

  const PetscScalar
          gridfactor = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dt_from_diffus = adaptTimeStepRatio
                     * 2 / ((gDmax + DEFAULT_ADDED_TO_GDMAX_ADAPT) * gridfactor);
  if (do_skip && (skipCountDown == 0)) {
    const PetscInt skip_max = static_cast<PetscInt>(config.get("skip_max"));
    const PetscScalar  conservativeFactor = 0.95;
    // typically "dt" in next line is from CFL for advection in temperature equation,
    //   but in fact it might be from other restrictions, e.g. CFL for mass continuity
    //   in basal sliding case, or max_dt
    skipCountDown = (PetscInt) floor(conservativeFactor * (dt / dt_from_diffus));
    skipCountDown = ( skipCountDown >  skip_max) ?  skip_max :  skipCountDown;
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

  bool do_mass_conserve = config.get_flag("do_mass_conserve"),
    do_energy = config.get_flag("do_energy");

  const PetscScalar timeToEnd = grid.time->end() - grid.time->current();
  if (dt_force > 0.0) {
    dt = dt_force; // override usual dt mechanism
    adaptReasonFlag = 'f';
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
  } else {
    dt = config.get("maximum_time_step_years", "years", "seconds");
    bool use_ssa_velocity = config.get_flag("use_ssa_velocity");

    adaptReasonFlag = 'm';

    if ((do_energy == PETSC_TRUE) && (doTemperatureCFL == PETSC_TRUE)) {
      // CFLmaxdt is set by computeMax3DVelocities() in call to velocity() iMvelocity.cc
      dt_from_cfl = CFLmaxdt;
      if (dt_from_cfl < dt) {
        dt = dt_from_cfl;
        adaptReasonFlag = 'c';
      }
    }
    if (btu) {
      PetscReal btu_dt;
      bool restrict;
      ierr = btu->max_timestep(grid.time->current(), btu_dt, restrict); CHKERRQ(ierr); // returns years
      if (restrict && btu_dt < dt) {
        dt = btu_dt;
        adaptReasonFlag = 'b';
      }
    }
    if (do_mass_conserve && use_ssa_velocity) {
      // CFLmaxdt2D is set by broadcastSSAVelocity()
      if (CFLmaxdt2D < dt) {
        dt = CFLmaxdt2D;
        adaptReasonFlag = 'u';
      }
    }
    if (do_mass_conserve) {
      // note: if do_skip then skipCountDown = floor(dt_from_cfl/dt_from_diffus)
      ierr = adaptTimeStepDiffusivity(); CHKERRQ(ierr); // might set adaptReasonFlag = 'd'
    }

    bool dteigencalving = config.get_flag("cfl_eigencalving");
    if (dteigencalving) {
      IceModelVec2V *ssa_velocity;
      ierr = stress_balance->get_2D_advective_velocity(ssa_velocity); CHKERRQ(ierr);
      ierr = stress_balance->compute_2D_principal_strain_rates(*ssa_velocity, vMask, strain_rates); CHKERRQ(ierr);
      ierr = dt_from_eigenCalving(); CHKERRQ(ierr);
      if (dt_from_eigencalving < dt) {
        dt = dt_from_eigencalving;
        adaptReasonFlag = 'k';
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


//! Because of the -skip mechanism it is still possible that we can have CFL violations: count them.
/*!
This applies to the horizontal part of the three-dimensional advection problem
solved by IceModel::ageStep() and the advection, ice-only part of the problem solved by
temperatureStep().  These methods use a fine vertical grid, and so we consider CFL
violations on that same fine grid. (FIXME: should we actually use the fine grid?)

Communication is needed to determine total CFL violation count over entire grid.
It is handled by temperatureAgeStep(), not here.
*/
PetscErrorCode IceModel::countCFLViolations(PetscScalar* CFLviol) {
  PetscErrorCode  ierr;

  const PetscScalar cflx = grid.dx / dt_TempAge,
                    cfly = grid.dy / dt_TempAge;

  PetscScalar *u, *v;
  IceModelVec3 *u3, *v3, *dummy;
  ierr = stress_balance->get_3d_velocity(u3, v3, dummy); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt  fks = grid.kBelowHeight(vH(i,j));

      ierr = u3->getInternalColumn(i,j,&u); CHKERRQ(ierr);
      ierr = v3->getInternalColumn(i,j,&v); CHKERRQ(ierr);

      // check horizontal CFL conditions at each point
      for (PetscInt k=0; k<=fks; k++) {
        if (PetscAbs(u[k]) > cflx)  *CFLviol += 1.0;
        if (PetscAbs(v[k]) > cfly)  *CFLviol += 1.0;
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = u3->end_access();  CHKERRQ(ierr);
  ierr = v3->end_access();  CHKERRQ(ierr);

  return 0;
}

