// Copyright (C) 2008-2014 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include "iceModelVec.hh"

#include "pism_options.hh"
#include <assert.h>

//! Finds PISM's input (-i or -boot_file) file using command-line options.
/*! This might be useful since coupling fields are usually in the file
  IceModel uses to initialize from.
*/
PetscErrorCode PISMComponent::find_pism_input(std::string &filename, //!< name of the file found
                                              bool &do_regrid, //!< specifies whether regridding is necessary
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
    do_regrid = true;
    start = 0;
  } else {
    do_regrid = false;
    start = last_record;
  }

  return 0;
}

/**
 * Regrid a variable by processing -regrid_file and -regrid_vars.
 *
 * @param[in] module_name Module name, used to annotate options when run with -help.

 * @param[out] variable pointer to an IceModelVec; @c variable has to
 *             have metadata set for this to work.
 *
 * @param[in] flag Regridding flag. If set to
 *            REGRID_WITHOUT_REGRID_VARS, regrid this variable by
 *            default, if =-regrid_vars= was not set. Otherwise a
 *            variable is only regridded if both =-regrid_file= and
 *            =-regrid_vars= are set *and* the name of the variable is
 *            found in the set of names given with =-regrid_vars=.
 *
 * @return 0 on success
 */
PetscErrorCode PISMComponent::regrid(std::string module_name, IceModelVec *variable,
                                     RegriddingFlag flag) {
  PetscErrorCode ierr;
  bool file_set, vars_set;
  std::set<std::string> vars;
  std::string file, title = module_name + std::string(" regridding options");

  assert(variable != NULL);

  ierr = PetscOptionsBegin(grid.com, "", title.c_str(), ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "regridding file name", file, file_set); CHKERRQ(ierr);
    ierr = PISMOptionsStringSet("-regrid_vars", "comma-separated list of regridding variables",
                                "", vars, vars_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (file_set == false)
    return 0;

  NCSpatialVariable &m = variable->metadata();

  if ((vars_set == true && set_contains(vars, m.get_string("short_name")) == true) ||
      (vars_set == false && flag == REGRID_WITHOUT_REGRID_VARS)) {
    ierr = verbPrintf(2, grid.com, "  regridding '%s' from file '%s' ...\n",
                      m.get_string("short_name").c_str(), file.c_str()); CHKERRQ(ierr);
    ierr = variable->regrid(file, CRITICAL); CHKERRQ(ierr);
  }

  return 0;
}
