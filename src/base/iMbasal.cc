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

//! Compute the coefficient for the basal velocity in SIA regions.
/*!
In SIA regions a basal sliding law of the form
  \f[ \mathbf{U}_b = (u_b,v_b) = - C \nabla h \f] 
is allowed.  Here \f$\mathbf{U}_b\f$ is the horizontal velocity of the base of the ice
(the "sliding velocity") and \f$h\f$ is the elevation of the ice surface.  This procedure 
returns the \em positive coefficient \f$C\f$ in this relationship.  This coefficient can
depend of the thickness, the basal temperature, and the horizontal location.

This procedure is virtual and can be replaced by any derived class.

The default version for IceModel here is location-independent 
pressure-melting-temperature-activated linear sliding.

THIS KIND OF SIA SLIDING LAW IS A BAD IDEA.  THAT'S WHY \f$\mu\f$ IS SET TO 
ZERO BY DEFAULT.                
 */
PetscScalar IceModel::basalVelocity(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {

  if (T + ice->beta_CC_grad * H > min_temperature_for_SIA_sliding) {
    return basal->velocity(mu, ice->rho * grav * H);
  } else {
    return 0;
  }
}


/*** for ice stream regions (MASK_DRAGGING): ***/
PetscScalar IceModel::basalDragx(PetscScalar **beta, PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(beta[i][j], tauc[i][j], u[i][j], v[i][j]);
}

PetscScalar IceModel::basalDragy(PetscScalar **beta, PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(beta[i][j], tauc[i][j], u[i][j], v[i][j]);
}


PetscErrorCode IceModel::initBasalTillModel() {
  PetscErrorCode ierr;
  
  if (createBasal_done == PETSC_FALSE) {
    if (doPlasticTill == PETSC_TRUE) {
      basal = new PlasticBasalType(plasticRegularization);
    } else {
      basal = new ViscousBasalType;
    }
    createBasal_done = PETSC_TRUE;
  }

  if (useSSAVelocity == PETSC_TRUE) {
    ierr = basal->printInfo(3,grid.com); CHKERRQ(ierr);
  }
  ierr = VecSet(vtauc, tauc_default_value); CHKERRQ(ierr);
  ierr = VecSet(vbeta, beta_default_drag_SSA); CHKERRQ(ierr);
  return 0;
}


//! Update the till yield stress for the plastic till model, based on pressure and stored till water.
/*! 
@cond CONTINUUM
We implement formula (2.4) in \lo\cite{SchoofStream}\elo.  That formula is
    \f[   \tau_c = \mu (\rho g H - p_w)\f]
We modify it by:

    (1) adding a small till cohesion \f$c_0\f$ (see \lo\cite{Paterson}\elo table 8.1);

    (2) replacing \f$p_w \to \lambda p_w\f$ where \f$\lambda =\f$ Hmelt / DEFAULT_MAX_HMELT;
       thus \f$0 \le \lambda \le 1\f$ always while \f$\lambda = 0\f$ when the bed is frozen; and

    (3) computing porewater pressure \f$p_w\f$ as a fixed fraction \f$\varphi\f$ of the overburden pressure \f$\rho g H\f$.

With these replacements our formula looks like
    \f[   \tau_c = c_0 + \mu \left(1 - \lambda \varphi\right) \rho g H \f]
Note also that \f$\mu = \tan(\theta)\f$ where \f$\theta\f$ is a ``friction angle''.  The parameters \f$c_0\f$, \f$\varphi\f$, \f$\theta\f$ can be set by options -till_cohesion, -till_pw_fraction, and -till_friction_angle, respectively.

@endcond
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
          const PetscScalar overburdenP = ice->rho * grav * H[i][j];
          const PetscScalar pwP = plastic_till_pw_fraction * overburdenP;
          // note Hmelt == 0 if frozen and  0 <= Hmelt <= Hmelt_max always
          //   so  0 <= lambda <= 1 
          const PetscScalar lambda = Hmelt[i][j] / Hmelt_max;
          const PetscScalar N = overburdenP - lambda * pwP;  // effective pressure on till
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

