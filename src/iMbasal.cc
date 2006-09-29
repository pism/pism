// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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
PetscScalar IceModel::basal(const PetscScalar H, const PetscScalar T, const PetscScalar mu) const {
  // This implements a pressure-melting-temperature activated linear sliding law.
  // returns *positive* coefficient C in the law:  U_b = <u_b,v_b> = - C grad h 
  if (T + ice.beta_CC_grad * H > DEFAULT_MIN_TEMP_FOR_SLIDING) {
    return mu * ice.rho * ice.grav * H;
  }
  return 0;
}

/*** for ice stream regions (MASK_DRAGGING): ***/
PetscScalar IceModel::basalDragx(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basalDrag(u[i][j], v[i][j]);
}

PetscScalar IceModel::basalDragy(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basalDrag(u[i][j], v[i][j]);
}

PetscScalar IceModel::basalDrag(const PetscScalar u, const PetscScalar v) const {
  /* This implements a linear sliding law.  A more elaborate relation may be
  * appropriate. */

  // return DEFAULT_BASAL_DRAG_COEFF_MACAYEAL;  // 2.0e9 Pa s m^-1
  return 4.0e9;  
}
