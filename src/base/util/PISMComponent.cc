// Copyright (C) 2008-2010 Ed Bueler and Constantine Khroulev
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

#include "PISMComponent.hh"
#include "PISMIO.hh"

//! Finds PISM's input (-i or -boot_file) file using command-line options.
/*! This might be useful since coupling fields are usually in the file
  IceModel uses to initialize from.
*/
PetscErrorCode PISMComponent::find_pism_input(string &filename, //!< name of the file found
					      LocalInterpCtx* &lic, //!< local interp. context
					      bool &regrid, //!< specifies whether regridding is necessary
					      int &start    //!< "start" to use when reading from filename
					      ) {
  PetscErrorCode ierr;
  PetscTruth i_set, boot_file_set;

  // read file names:
  char i_file[PETSC_MAX_PATH_LEN], boot_file_file[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-i", i_file, 
			       PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_file", boot_file_file, 
			       PETSC_MAX_PATH_LEN, &boot_file_set); CHKERRQ(ierr);
  if (i_set) {
    if (boot_file_set) {
      ierr = PetscPrintf(grid.com,
	"PISMClimateCoupler ERROR: both '-i' and '-boot_file' are used. Exiting...\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    filename = i_file;
  }
  else if (boot_file_set) {
    filename = boot_file_file;
  }

  // filename now contains name of PISM input (or bootstrapping) file; now check
  // it is really there; if so, read the dimensions of computational grid so
  // that we can set up a LocalInterpCtx for actual reading of climate data
  PISMIO nc(&grid);
  int last_record;
  grid_info gi;
  ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
  last_record -= 1;
  ierr = nc.close(); CHKERRQ(ierr);

  if (boot_file_set) {
    // *caller* of find_pism_input() is in charge of destroying
    lic = new LocalInterpCtx(gi, NULL, NULL, grid); // 2D only
    regrid = true;
    start = 0;
  } else {
    lic = NULL;
    regrid = false;
    start = last_record;
  }

  return 0;
}
