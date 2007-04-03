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

/*
* The basal sliding components. These could be wrapped in a BasalType class ala IceType.
* That would allow easy plugability of different sliding relations. It should be
* very easy to implement when an alternative sliding law is included.
*/

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

PetscInt IceModel::initBasalFields() {
  PetscErrorCode ierr;
  ierr = VecSet(vtauc, DEFAULT_TAUC); CHKERRQ(ierr);
  ierr = VecSet(vbeta, DEFAULT_BASAL_DRAG_COEFF_MACAYEAL); CHKERRQ(ierr);
  return 0;
}
