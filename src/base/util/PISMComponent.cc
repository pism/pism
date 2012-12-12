// Copyright (C) 2008-2012 Ed Bueler and Constantine Khroulev
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
#include "PIO.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "NCVariable.hh"

//! Finds PISM's input (-i or -boot_file) file using command-line options.
/*! This might be useful since coupling fields are usually in the file
  IceModel uses to initialize from.
*/
PetscErrorCode PISMComponent::find_pism_input(string &filename, //!< name of the file found
					      bool &regrid, //!< specifies whether regridding is necessary
					      int &start    //!< "start" to use when reading from filename
					      ) {
  PetscErrorCode ierr;
  PetscBool i_set, boot_file_set;

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
      PISMEnd();
    }
    filename = i_file;
  }
  else if (boot_file_set) {
    filename = boot_file_file;
  }

  PIO nc(grid, "netcdf3");      // OK to use netcdf3
  unsigned int last_record;
  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
  ierr = nc.inq_nrecords(last_record); CHKERRQ(ierr);
  last_record -= 1;
  ierr = nc.close(); CHKERRQ(ierr);

  if (boot_file_set) {
    regrid = true;
    start = 0;
  } else {
    regrid = false;
    start = last_record;
  }

  return 0;
}
