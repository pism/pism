// Copyright (C) 2009 Constantine Khroulev
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

//! Initializes the code writing scalar time-series.
PetscErrorCode IceModel::init_scalar_timeseries() {
  PetscErrorCode ierr;
  PetscTruth ts_file = PETSC_FALSE, ts_times = PETSC_FALSE;
  char tmp[TEMPORARY_STRING_LENGTH] = "\0";

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_filename", tmp,
			       PETSC_MAX_PATH_LEN, &ts_file); CHKERRQ(ierr);
  scalar_ts_filename = tmp;

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_times", tmp,
			       TEMPORARY_STRING_LENGTH, &ts_times); CHKERRQ(ierr);

  if (ts_file ^ ts_times) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: you need to specity both -ts_filename and -ts_times to save"
		       "diagnostic time-seties.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // If neither -ts_filename nor -ts_times is set, we're done.
  if (!ts_file && !ts_times) {
    save_scalar_ts = false;
    return 0;
  }
  
  save_scalar_ts = true;

  ierr = parse_times(grid.com, tmp, scalar_ts_times);
  if (ierr != 0) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: parsing the -ts_times argument failed.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = verbPrintf(2, grid.com, "saving scalar time-series to '%s'; ",
		    scalar_ts_filename.c_str()); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "times requested: %s\n", tmp); CHKERRQ(ierr);

  current_scalar_ts = 0;
  return 0;
}


PetscErrorCode IceModel::write_scalar_timeseries() {
  PetscErrorCode ierr;

  // return if no time-series requested
  if (!save_scalar_ts) return 0;

  // return if wrote all the records already
  if (current_scalar_ts == scalar_ts_times.size())
    return 0;

  // return if did not yet reach the time we need to save at
  if (scalar_ts_times[current_scalar_ts] > grid.year)
    return 0;

  // compute_scalar_timeseries();

  while ((scalar_ts_times[current_scalar_ts] <= grid.year) &&
	 (current_scalar_ts < scalar_ts_times.size())) {
    
    // write_diagnostics();

    ierr = PetscPrintf(grid.com, "Would save at t = %f\n",
		       scalar_ts_times[current_scalar_ts]); CHKERRQ(ierr);

    current_scalar_ts++;
  }

  return 0;
}
