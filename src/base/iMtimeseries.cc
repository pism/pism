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
  bool save_at_equal_intervals = false;
  double a, delta, b;
  PetscInt N;

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_filename", tmp,
			       PETSC_MAX_PATH_LEN, &ts_file); CHKERRQ(ierr);
  scalar_ts_filename = tmp;

  ierr = PetscOptionsGetString(PETSC_NULL, "-ts_times", tmp,
			       TEMPORARY_STRING_LENGTH, &ts_times); CHKERRQ(ierr);

  if (ts_file && !ts_times) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: -ts_filename is set, but -ts_times is not.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  if (ts_times && !ts_file) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: -ts_times is set, but -ts_filename is not.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // If neither -ts_filename nor -ts_times is set, we're done.
  if (!ts_file && !ts_times) {
    save_scalar_ts = false;
    return 0;
  } else
    save_scalar_ts = true;

  if (strchr(tmp, ':')) {	// if a string contains a colon...

    ierr = parse_range(grid.com, tmp, &a, &delta, &b);

    if (ierr != 0) {
      ierr = PetscPrintf(grid.com,
         "PISM ERROR: Parsing the -ts_times argument (%s) failed.\n",
	 tmp); CHKERRQ(ierr);
      PetscEnd();
    }

    if (a >= b) {
      ierr = PetscPrintf(grid.com,
         "PISM ERROR: Error in the -ts_times argument: a >= b in the range specification '%s'.\n",
	 tmp); CHKERRQ(ierr);
      PetscEnd();
    }

    if (delta <= 0) {
      ierr = PetscPrintf(grid.com,
	 "PISM ERROR: Error in the -ts_times argument: dt <= 0 in the range specification '%s'.\n",
	 tmp); CHKERRQ(ierr);
      PetscEnd();
    }

    // Compute the number of snapshots and the times to save after
    N = floor((b - a)/delta) + 1; // number of snapshots
    scalar_ts_times.resize(N);

    for (int j = 0; j < N; ++j)
      scalar_ts_times[j] = a + delta*j;

    save_at_equal_intervals = true;
  } else {			// no colon; it must be a list of numbers
    N = (int)config.get("max_number_of_snapshots");
    scalar_ts_times.resize(N);

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-ts_times", &scalar_ts_times[0], &N, PETSC_NULL);
    scalar_ts_times.resize(N);

    sort(scalar_ts_times.begin(), scalar_ts_times.end());
  }

  ierr = verbPrintf(2, grid.com, "saving snapshots to '%s'; ",
		    snapshots_filename.c_str()); CHKERRQ(ierr);


  return 0;
}

PetscErrorCode IceModel::write_scalar_timeseries() {

  return 0;
}
