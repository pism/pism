// Copyright (C) 2010 Constantine Khroulev and Ed Bueler
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

#include "PISMStressBalance.hh"

//! \brief Performs the shallow stress balance computation.
PetscErrorCode PISMStressBalance::update(bool fast) {
  PetscErrorCode ierr;
  IceModelVec2V *velocity_2d;
  IceModelVec3  *u, *v;

  ierr = stress_balance->update(fast); CHKERRQ(ierr);

  ierr = stress_balance->get_advective_2d_velocity(velocity_2d); CHKERRQ(ierr); 

  ierr = modifier->update(velocity_2d, fast); CHKERRQ(ierr);

  if (!fast) {
    ierr = modifier->get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);

    ierr = compute_vertical_velocity(u, v); CHKERRQ(ierr); 
  }

  return 0;
}
