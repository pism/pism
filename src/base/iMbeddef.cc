// Copyright (C) 2004-2011 Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"

PetscErrorCode IceModel::allocate_bed_deformation() {
  PetscErrorCode ierr;
  string model = config.get_string("bed_deformation_model");
  set<string> choices;

  ierr = check_old_option_and_stop(grid.com, "-bed_def_iso", "-bed_def"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-bed_def_lc",  "-bed_def"); CHKERRQ(ierr);
  
  choices.insert("none");
  choices.insert("iso");
#if (PISM_HAVE_FFTW==1)
  choices.insert("lc");
#endif

  ierr = PetscOptionsBegin(grid.com, "", "Bed deformation model", ""); CHKERRQ(ierr);
  {
    bool dummy;
    ierr = PISMOptionsList(grid.com, "-bed_def", "Specifies a bed deformation model.",
			 choices, model, model, dummy); CHKERRQ(ierr);

  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (model == "none")
    return 0;

  if ((model == "iso") && (beddef == NULL)) {
    beddef = new PBPointwiseIsostasy(grid, config);
    return 0;
  }

#if (PISM_HAVE_FFTW==1)
  if ((model == "lc") && (beddef == NULL)) {
    beddef = new PBLingleClark(grid, config);
    return 0;
  }
#endif

  return 0;
}

PetscErrorCode IceModel::bed_def_step(bool &bed_changed) {
  PetscErrorCode  ierr;

  if (beddef == NULL) SETERRQ(1, "beddef == NULL");

  double update_interval = config.get("bed_def_interval_years", "years", "seconds");

  // This is a front end to the bed deformation update system.  It updates
  // no more often than update_interval.

  // If the bed elevations are not expired, exit cleanly.
  const PetscScalar dt_beddef = grid.time->current() - last_bed_def_update;
  if (dt_beddef >= update_interval) {

    ierr = beddef->update(grid.time->current(), dt); CHKERRQ(ierr);

    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

    last_bed_def_update = grid.time->current();
    // Mark bed topography as "modified":
    vbed.inc_state_counter();
    // SIAFD::update() uses this to decide if we need to re-compute the
    // smoothed bed.
    bed_changed = true;
  } else
    bed_changed = false;

  return 0;
}

