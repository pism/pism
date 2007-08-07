// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"

/*** for SIA regions (MASK_SHEET): ***/
PetscScalar IceModel::basalVelocity(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {

  //PetscErrorCode  ierr = PetscPrintf(grid.com, 
  //        "   [IceModel::basal called with:   x=%f, y=%f, H=%f, T=%f, alpha=%f]\n",
  //        x,y,H,T,alpha);  CHKERRQ(ierr);

  // This implements location-independent pressure-melting-temperature-activated
  // linear sliding law.  Returns *positive* coefficient C in the law
  //                U_b = <u_b,v_b> = - C grad h 
  if (T + ice.beta_CC_grad * H > DEFAULT_MIN_TEMP_FOR_SLIDING) {
    return basal->velocity(mu, ice.rho * grav * H);
  }
  return 0;
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
  if (useMacayealVelocity == PETSC_TRUE) {
    ierr = basal->printInfo(3,grid.com); CHKERRQ(ierr);
  }
  ierr = VecSet(vtauc, DEFAULT_TAUC); CHKERRQ(ierr);
  ierr = VecSet(vbeta, DEFAULT_BASAL_DRAG_COEFF_MACAYEAL); CHKERRQ(ierr);
  createBasal_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceModel::updateYieldStressFromHmelt() {
  PetscErrorCode  ierr;
  // only makes sense when doPlasticTill == TRUE
  // we implement formula (2.4) in C. Schoof 2006 "A variational approach to ice
  // stream flow", J. Fluid Mech. vol 556 pp 227--251:
  //  (2.4)   \tau_c = \mu (\rho g H - p_w)
  // we modify it by:
  //   1. adding a small till cohesion (see Paterson 3rd ed table 8.1)
  //   2. replacing   p_w --> \lambda p_w   where \lambda = Hmelt / DEFAULT_MAX_HMELT;
  //      thus 0 <= \lambda <= 1 and \lambda = 0 when bed is frozen 
  //   3. computing a porewater pressure p_w which is the max of 0.95 * overburden
  //      and the porewater pressure computed by formula (4) in 
  //      C. Ritz et al 2001 J. G. R. vol 106 no D23 pp 31943--31964;
  //      the modification of this porewater pressure as in Lingle&Brown 1987 is not 
  //      implementable because the "elevation of the bed at the grounding line"
  //      is at an unknowable location (we are not doing a flow line model!)

//  const PetscScalar plastic_till_c_0 = 20.0e3;  // Pa; 20kPa = 0.2 bar; cohesion of till
  const PetscScalar plastic_till_c_0 = 5.0e3;
//  const PetscScalar plastic_till_mu = 0.466307658156;  // = tan(25^o); till friction angle
  const PetscScalar plastic_till_mu = 0.2125565616700221;  // = tan(12^o); till friction angle
    
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
//          const PetscScalar drivingP = - ocean.rho * grav * bed[i][j];
//          const PetscScalar pw = PetscMax(porewater_gamma * overburdenP, drivingP);
//          const PetscScalar pw = porewater_gamma * overburdenP;
        const PetscScalar bedfrac = PetscMax(-bed[i][j],0.0) / 1000.0;
        const PetscScalar pw = (0.85 + 0.1 * PetscMin(bedfrac,1.0)) * overburdenP;
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
  // communicate mask
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}

