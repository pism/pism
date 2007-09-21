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

The default version for IceModel here is location-independent pressure-melting-temperature-activated 
linear sliding.
  //                
 */
PetscScalar IceModel::basalVelocity(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {

  if (T + ice.beta_CC_grad * H > DEFAULT_MIN_TEMP_FOR_SLIDING) {
    return basal->velocity(mu, ice.rho * grav * H);
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
      basal = new PlasticBasalType;
    } else {
      basal = new ViscousBasalType;
    }
  }
  if (useSSAVelocity == PETSC_TRUE) {
    ierr = basal->printInfo(3,grid.com); CHKERRQ(ierr);
  }
  ierr = VecSet(vtauc, DEFAULT_TAUC); CHKERRQ(ierr);
  ierr = VecSet(vbeta, DEFAULT_BASAL_DRAG_COEFF_SSA); CHKERRQ(ierr);
  createBasal_done = PETSC_TRUE;
  return 0;
}


//! Update the till yield stress for the plastic till model, based on pressure and stored till water.
/*! 
@cond CONTINUUM
We implement formula (2.4) in C. Schoof 2006 "A variational approach to ice stream flow", J. Fluid Mech. vol 556 pp 227--251.  That formula is
    \f[   \tau_c = \mu (\rho g H - p_w)\f]
We modify it by:

    (1) adding a small till cohesion \f$c_0\f$ (see Paterson 3rd ed table 8.1);

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

  PetscScalar **mask, **tauc, **H, **Hmelt, **bed; 
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        tauc[i][j] = 0.0;  
      } else { // grounded
        mask[i][j] = MASK_DRAGGING;  // in Schoof model, everything is dragging, so force this
        const PetscScalar overburdenP = ice.rho * grav * H[i][j];
#if 0
//  const PetscScalar plastic_till_c_0 = 20.0e3;  // Pa; 20kPa = 0.2 bar; cohesion of till
//  plastic_till_c_0 = 5.0e3;
//  const PetscScalar plastic_till_mu = 0.466307658156;  // = tan(25^o); till friction angle
//  plastic_till_mu = 0.2125565616700221;  // = tan(12^o); till friction angle
//          const PetscScalar drivingP = - ocean.rho * grav * bed[i][j];
//          const PetscScalar pw = PetscMax(porewater_gamma * overburdenP, drivingP);
//          const PetscScalar pw = porewater_gamma * overburdenP;
        const PetscScalar bedfrac = PetscMax(-bed[i][j],0.0) / 1000.0;
        const PetscScalar pw = (0.85 + 0.1 * PetscMin(bedfrac,1.0)) * overburdenP;
#else
        const PetscScalar pw = plastic_till_pw_fraction * overburdenP;
#endif
        const PetscScalar lambda = Hmelt[i][j] / DEFAULT_MAX_HMELT;  // note Hmelt[i][j]=0 if frozen
        tauc[i][j] = plastic_till_c_0 + plastic_till_mu * (overburdenP - lambda * pw);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  // communicate possibly updated mask
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}

