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

PetscErrorCode SSB_Modifier::init(PISMvars &vars) {
  PetscErrorCode ierr;

  ierr =     u.create(grid, "uvel", true); CHKERRQ(ierr);
  ierr =     u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
			  "m s-1", "land_ice_x_velocity"); CHKERRQ(ierr);
  ierr =     u.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  u.write_in_glaciological_units = true;

  ierr =     v.create(grid, "vvel", true); CHKERRQ(ierr);
  ierr =     v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
			  "m s-1", "land_ice_y_velocity"); CHKERRQ(ierr);
  ierr =     v.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  v.write_in_glaciological_units = true;

  ierr = Sigma.create(grid, "strainheat", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr = Sigma.set_attrs("internal",
                          "rate of strain heating in ice (dissipation heating)",
	        	  "W m-3", ""); CHKERRQ(ierr);
  ierr = Sigma.set_glaciological_units("mW m-3"); CHKERRQ(ierr);

  ierr = diffusive_flux.create(grid, "diffusive_flux", true, 1); CHKERRQ(ierr);
  ierr = diffusive_flux.set_attrs("internal", 
                                  "diffusive (SIA) flux components on the staggered grid",
                                  "", ""); CHKERRQ(ierr);

  return 0; 
}

PetscErrorCode SSB_Modifier::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;

  ierr =     u.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr =     v.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = Sigma.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);

  return 0;
}



//! \brief Distribute the input velocity throughout the column.
PetscErrorCode SSBM_Trivial::update(IceModelVec2V &vel_input, bool fast) {
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
PetscErrorCode SSBM_Trivial::compute_sigma(IceModelVec2S *D2_input, IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  
  ierr = enthalpy->end_access(); CHKERRQ(ierr);

  return 0;
}
