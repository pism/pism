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

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"


//! Compute the coefficient of surface gradient, for basal sliding velocity as a function of driving stress in SIA regions.
/*!
THIS KIND OF SIA SLIDING LAW IS A BAD IDEA.  THAT'S WHY \f$\mu\f$ IS SET TO 
ZERO BY DEFAULT.                

In SIA regions a basal sliding law of the form
  \f[ \mathbf{U}_b = (u_b,v_b) = - C \nabla h \f] 
is allowed.  Here \f$\mathbf{U}_b\f$ is the horizontal velocity of the base of the ice
(the "sliding velocity") and \f$h\f$ is the elevation of the ice surface.  This procedure 
returns the \em positive coefficient \f$C\f$ in this relationship.  This coefficient can
depend of the thickness, the basal temperature, and the horizontal location.

This procedure is virtual and can be replaced by any derived class.

The default version for IceModel here is location-independent 
pressure-melting-temperature-activated linear sliding.  Here we pass
\f$\mu\f$, which can be set by option \c -mu_sliding, and the pressure 
at the base to BasalTypeSIA::velocity().

The returned coefficient is used in basalSlidingHeatingSIA().
 */
PetscScalar IceModel::basalVelocity(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {

  if (T + ice->beta_CC_grad * H > min_temperature_for_SIA_sliding) {
    return basalSIA->velocity(mu, ice->rho * grav * H);
  } else {
    return 0;
  }
}


/*** for ice stream regions (MASK_DRAGGING): ***/
PetscScalar IceModel::basalDragx(PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], u[i][j], v[i][j]);
}

PetscScalar IceModel::basalDragy(PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], u[i][j], v[i][j]);
}


//! Initialize the pseudo-plastic till mechanical model.
/*! 
See PlasticBasalType and updateYieldStressFromHmelt() and getEffectivePressureOnTill()
for model equations.  See also invertVelocitiesFromNetCDF() for one way to 
get a map of till friction angle.

Also initializes the SIA-type sliding law, but use of that model is not recommended
and is turned off by default.
 */
PetscErrorCode IceModel::initBasalTillModel() {
  PetscErrorCode ierr;
  
  if (createBasal_done == PETSC_FALSE) {
    basal = new PlasticBasalType(plasticRegularization, doPseudoPlasticTill, 
                                 pseudo_plastic_q, pseudo_plastic_uthreshold);
    basalSIA = new BasalTypeSIA();  // initialize it; USE NOT RECOMMENDED!
    createBasal_done = PETSC_TRUE;
  }

  if (useSSAVelocity == PETSC_TRUE) {
    ierr = basal->printInfo(3,grid.com); CHKERRQ(ierr);
  }
  ierr = VecSet(vtauc, tauc_default_value); CHKERRQ(ierr);
  // since vtillphi is part of model state it should not be set to default here, but 
  //   rather as part of initialization/bootstrapping
  return 0;
}


//! Compute effective pressure on till using effective thickness of stored till water.
/*!
Uses ice thickness to compute overburden pressure.  Pore water pressure is assumed
to be a fixed fraction of the overburden pressure.

Note \c melt_thk should be zero at points where base of ice is frozen.

Also we always want \f$0 \le\f$ \c melt_thk \f$\le\f$ \c Hmelt_max 
so \f$0 \le\f$ \c lambda \f$\le 1\f$ inside this routine.
 */
PetscScalar IceModel::getEffectivePressureOnTill(
               const PetscScalar thk, const PetscScalar melt_thk) {
  const PetscScalar
     overburdenP = ice->rho * grav * thk,
     pwP = plastic_till_pw_fraction * overburdenP,
     lambda = melt_thk / Hmelt_max;
  return overburdenP - lambda * pwP;  
}


//! Update the till yield stress for the pseudo-plastic till model.
/*!
Expanded brief description: Update the till yield stress \e and the mask, for 
the pseudo-plastic till model, based on pressure and stored till water.

This procedure also modifies the mask.  In particular, it has the side effect
of marking all grounded points as MASK_DRAGGING.  (FIXME:  This aspect should be
refactored.  Unnecessary communication can probably be avoided.)

We implement formula (2.4) in \lo\cite{SchoofStream}\elo.  That formula is
    \f[   \tau_c = \mu (\rho g H - p_w)\f]
We modify it by:
    - adding a small till cohesion \f$c_0\f$ (see \lo\cite{Paterson}\elo table 8.1);
    - replacing \f$p_w \to \lambda p_w\f$ where \f$\lambda =\f$ 
      Hmelt/DEFAULT_MAX_HMELT; thus \f$0 \le \lambda \le 1\f$ always while 
      \f$\lambda = 0\f$ when the bed is frozen; and
    - computing porewater pressure \f$p_w\f$ as a fixed fraction \f$\varphi\f$ 
      of the overburden pressure \f$\rho g H\f$.
The effective pressure \f$\rho g H - p_w\f$ is actually computed by 
getEffectivePressureOnTill().

With these replacements our formula looks like
    \f[   \tau_c = c_0 + \mu \left(1 - \lambda \varphi\right) \rho g H \f]
Note also that \f$\mu = \tan(\theta)\f$ where \f$\theta\f$ is a "friction angle."
The parameters \f$c_0\f$, \f$\varphi\f$, \f$\theta\f$ can be set by options 
\c -till_cohesion, \c -till_pw_fraction, and \c -till_friction_angle, respectively.
 */
PetscErrorCode IceModel::updateYieldStressFromHmelt() {
  PetscErrorCode  ierr;
  //      (compare the porewater pressure computed by formula (4) in 
  //      C. Ritz et al 2001 J. G. R. vol 106 no D23 pp 31943--31964;
  //      the modification of this porewater pressure as in Lingle&Brown 1987 is not 
  //      implementable because the "elevation of the bed at the grounding line"
  //      is at an unknowable location as we are not doing a flow line model!)

  // only makes sense when doPlasticTill == TRUE
  if (doPlasticTill == PETSC_FALSE) {
    SETERRQ(1,"doPlasticTill == PETSC_FALSE but updateYieldStressFromHmelt() called");
  }

  if (holdTillYieldStress == PETSC_TRUE) {  // don't modify tauc; use stored
    PetscScalar **mask; 
    ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (modMask(mask[i][j]) != MASK_FLOATING) {
          mask[i][j] = MASK_DRAGGING;  // in Schoof model, everything grounded is dragging,
                                       // so force this
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  } else { // usual case: use Hmelt to determine tauc
    PetscScalar **mask, **tauc, **H, **Hmelt, **bed, **tillphi; 
    ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vtillphi, &tillphi); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          tauc[i][j] = 0.0;  
        } else if (H[i][j] == 0.0) {
          tauc[i][j] = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar if no ice
          mask[i][j] = MASK_DRAGGING;  // mark it this way anyway
        } else { // grounded and there is some ice
          mask[i][j] = MASK_DRAGGING;  // in Schoof model, everything is dragging, so force this
          const PetscScalar N = getEffectivePressureOnTill(H[i][j], Hmelt[i][j]);
          if (useConstantTillPhi == PETSC_TRUE) {
            tauc[i][j] = plastic_till_c_0 + plastic_till_mu * N;
          } else {
            const PetscScalar mymu = tan((pi/180.0) * tillphi[i][j]);
            tauc[i][j] = plastic_till_c_0 + mymu * N;
          }
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vtillphi, &tillphi); CHKERRQ(ierr);
  }

  // communicate possibly updated mask; tauc does not need communication
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}


//! Apply explicit time step for pure diffusion to basal layer of melt water.
/*!
See preprint \lo\cite{BBssasliding}\elo.

Uses vWork2d[0] to temporarily store new values for Hmelt.
 */
PetscErrorCode IceModel::diffuseHmelt() {
  PetscErrorCode  ierr;
  
  // diffusion constant K in u_t = K \nabla^2 u is chosen so that fundmental
  //   solution has standard deviation \sigma = 20 km at time t = 1000 yrs;
  //   2 \sigma^2 = 4 K t
  const PetscScalar K = 2.0e4 * 2.0e4 / (2.0 * 1000.0 * secpera),
                    Rx = K * dtTempAge / (grid.dx * grid.dx),
                    Ry = K * dtTempAge / (grid.dy * grid.dy);

  // NOTE: restriction that
  //    1 - 2 R_x - 2 R_y \ge 0
  // is a maximum principle restriction; therefore new Hmelt will be between
  // zero and Hmelt_max if old Hmelt has that property
  const PetscScalar oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  if (oneM4R <= 0.0) {
    SETERRQ(1,
       "diffuseHmelt() has 1 - 2Rx - 2Ry <= 0 so explicit method for diffusion unstable\n"
       "  (timestep restriction believed so rare that is not part of adaptive scheme)");
  }

  // communicate ghosted values so neighbors are valid
  ierr = DALocalToLocalEnd(grid.da2, vHmelt, INSERT_VALUES, vHmelt); CHKERRQ(ierr);

  PetscScalar **Hmelt, **Hmeltnew; 
  ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &Hmeltnew); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Hmeltnew[i][j] = oneM4R * Hmelt[i][j]
                       + Rx * (Hmelt[i+1][j] + Hmelt[i-1][j])
                       + Ry * (Hmelt[i][j+1] + Hmelt[i][j-1]);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &Hmeltnew); CHKERRQ(ierr);

  // finally copy new into vHmelt (and communicate ghosted values at same time)
  ierr = DALocalToLocalBegin(grid.da2, vWork2d[0], INSERT_VALUES, vHmelt); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vWork2d[0], INSERT_VALUES, vHmelt); CHKERRQ(ierr);

  return 0;
}

