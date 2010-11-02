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

#include "SSB_Modifier.hh"

//! \brief Distribute the input velocity throughout the column.
PetscErrorCode SSB_Trivial::update(IceModelVec2V &vel_input, bool fast) {
  PetscErrorCode ierr;

  ierr = u.begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = vel_input.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u.setColumn(i,j, vel_input(i,j).u); CHKERRQ(ierr);
      ierr = v.setColumn(i,j, vel_input(i,j).v); CHKERRQ(ierr);
    }
  }

  ierr = vel_input.end_access(); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);  

  return 0;
}


//! \brief Compute the volumetric strain heating.
/*!
 * Uses 
 * - delta on the staggered grid, which should be initialized by the update(true) call.
 * - enthalpy
 * - surface gradient on the staggered grid
 * - ice thickness relative to the smoothed bed
 */
PetscErrorCode SSB_Modifier::compute_sigma(IceModelVec2S &D2_input, IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  
  ierr = enthalpy->end_access(); CHKERRQ(ierr);

  return 0;
}
