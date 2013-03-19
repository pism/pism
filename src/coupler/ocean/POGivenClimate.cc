// Copyright (C) 2011, 2012, 2013 Constantine Khroulev
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

#include "POGivenClimate.hh"
#include "IceGrid.hh"

POGiven::POGiven(IceGrid &g, const NCConfigVariable &conf)
  : PGivenClimate<POModifier,PISMOceanModel>(g, conf, NULL)
{
  PetscErrorCode ierr = allocate_POGiven(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

POGiven::~POGiven() {
  // empty
}

PetscErrorCode POGiven::allocate_POGiven() {
  PetscErrorCode ierr;
  option_prefix   = "-ocean_given";

  // will be de-allocated by the parent's destructor
  shelfbtemp     = new IceModelVec2T;
  shelfbmassflux = new IceModelVec2T;

  m_fields["shelfbtemp"]     = shelfbtemp;
  m_fields["shelfbmassflux"] = shelfbmassflux;

  ierr = process_options(); CHKERRQ(ierr);

  map<string, string> standard_names;
  ierr = set_vec_parameters(standard_names); CHKERRQ(ierr);

  ierr = shelfbtemp->create(grid, "shelfbtemp", false); CHKERRQ(ierr);
  ierr = shelfbmassflux->create(grid, "shelfbmassflux", false); CHKERRQ(ierr);

  ierr = shelfbtemp->set_attrs("climate_forcing",
                        "absolute temperature at ice shelf base",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = shelfbmassflux->set_attrs("climate_forcing",
			     "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
			     "m s-1", ""); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGiven::init(PISMVars &) {
  PetscErrorCode ierr;

  t = dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean model reading base of the shelf temperature\n"
                    "  and sub-shelf mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = shelfbtemp->init(filename); CHKERRQ(ierr);
  ierr = shelfbmassflux->init(filename); CHKERRQ(ierr);

  // read time-independent data right away:
  if (shelfbtemp->get_n_records() == 1 && shelfbmassflux->get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  return 0;
}

PetscErrorCode POGiven::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = shelfbmassflux->at_time(t, bc_period, bc_reference_time); CHKERRQ(ierr);
  ierr = shelfbtemp->at_time(t, bc_period, bc_reference_time); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGiven::sea_level_elevation(PetscReal &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POGiven::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = shelfbtemp->copy_to(result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode POGiven::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = shelfbmassflux->copy_to(result); CHKERRQ(ierr);
  return 0;
}

