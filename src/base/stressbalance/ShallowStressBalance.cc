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
PetscErrorCode ShallowStressBalance::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;

  ierr = velocity.create(grid, "bar", true); CHKERRQ(ierr); // components are ubar and vbar
  ierr = velocity.set_attrs("model_state", "thickness-advective ice velocity", "m s-1", ""); CHKERRQ(ierr);
  ierr = velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  velocity.write_in_glaciological_units = true;

  ierr = basal_frictional_heating.create(grid, "bfrict", false); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_attrs("diagnostic", "", "W m-2", ""); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_glaciological_units("mW m-2"); CHKERRQ(ierr);
  basal_frictional_heating.write_in_glaciological_units = true;

  ierr = D2.create(grid, "D2", false); CHKERRQ(ierr);
  ierr = D2.set_attrs("internal",
                      "(partial) square of the Frobenius norm of D_{ij}, the combined strain rates",
                      "", ""); CHKERRQ(ierr);

  return 0;
}

//! \brief Update the trivial shallow stress balance object.
PetscErrorCode SSB_Trivial::update(bool fast) {
  PetscErrorCode ierr;
  if (fast) return 0;

  ierr = velocity.set(0.0); CHKERRQ(ierr);
  max_u = max_v = 0.0;

  ierr = basal_frictional_heating.set(0.0); CHKERRQ(ierr);

  ierr = D2.set(0.0); CHKERRQ(ierr);
  
  return 0;
}

