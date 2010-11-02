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

#include "ShallowStressBalance.hh"

//! \brief Initialize a shallow stress balance object.
PetscErrorCode ShallowStressBalance::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = velocity.create(grid, "bar", true); CHKERRQ(ierr); // ubar and vbar
  ierr = velocity.set_attrs("model_state", "thickness-advective ice velocity", "m s-1", ""); CHKERRQ(ierr);
  velocity.time_independent = false;
  ierr = velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  velocity.write_in_glaciological_units = true;
  ierr = variables.add(velocity); CHKERRQ(ierr);

  ierr = basal_frictional_heating.create(grid, "bfrict", false); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_attrs("diagnostic", "", "W m-2", ""); CHKERRQ(ierr);
  basal_frictional_heating.time_independent = false;
  ierr = basal_frictional_heating.set_glaciological_units("mW m-2"); CHKERRQ(ierr);
  basal_frictional_heating.write_in_glaciological_units = true;
  ierr = variables.add(basal_frictional_heating); CHKERRQ(ierr);

  max_u = max_v = 0;

  return 0;
}

//! \brief Specify velocity boundary conditions.
PetscErrorCode ShallowStressBalance::set_boundary_conditions(IceModelVec2Mask &locations,
                                                             IceModelVec2V &bc) {
  set_bc = true;
  bc_locations = &locations;
  vel_bc = &bc;
  return 0;
}

//! \brief Get a pointer to the thickness-advective 2D velocity field.
PetscErrorCode ShallowStressBalance::get_advective_2d_velocity(IceModelVec2V* &result) {
  result = &velocity;
  return 0;
}

//! \brief Get the maximum x- and y-components of the 2D velocity.
PetscErrorCode ShallowStressBalance::get_max_2d_velocity(PetscReal &u_max, PetscReal &v_max) {
  u_max = max_u;
  v_max = max_v;
  return 0;
}

//! \brief Get the basal frictional heating field.
PetscErrorCode ShallowStressBalance::get_basal_frictional_heating(IceModelVec2S* &result) {
  result = &basal_frictional_heating;
  return 0;
}
