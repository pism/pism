// Copyright (C) 2008-2015 Ed Bueler and Constantine Khroulev
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

#include "error_handling.hh"

namespace pism {

const IceGrid& Component::get_grid() const {
  return m_grid;
}

//! Finds PISM's input (-i or -boot_file) file using command-line options.
/*! This might be useful since coupling fields are usually in the file
  IceModel uses to initialize from.
*/
void Component::find_pism_input(std::string &filename, bool &do_regrid, int &start) {
  PetscErrorCode ierr;
  PetscBool i_set, boot_file_set;

  // read file names:
  char i_file[PETSC_MAX_PATH_LEN], boot_file_file[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(NULL, "-i", i_file, 
                               PETSC_MAX_PATH_LEN, &i_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetString");
  ierr = PetscOptionsGetString(NULL, "-boot_file", boot_file_file, 
                               PETSC_MAX_PATH_LEN, &boot_file_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetString");
  if (i_set) {
    if (boot_file_set) {
      throw RuntimeError("both '-i' and '-boot_file' are used.");
    }
    filename = i_file;
  }
  else if (boot_file_set) {
    filename = boot_file_file;
  }

  PIO nc(m_grid, "netcdf3");      // OK to use netcdf3
  unsigned int last_record;
  nc.open(filename, PISM_READONLY);
  last_record = nc.inq_nrecords() - 1;
  nc.close();

  if (boot_file_set) {
    do_regrid = true;
    start = 0;
  } else {
    do_regrid = false;
    start = last_record;
  }
}

/**
 * Regrid a variable by processing -regrid_file and -regrid_vars.
 *
 * @param[in] module_name Module name, used to annotate options when run with -help.
 *
 * @param[out] variable pointer to an IceModelVec; @c variable has to
 *             have metadata set for this to work.
 *
 * @param[in] flag Regridding flag. If set to
 *            REGRID_WITHOUT_REGRID_VARS, regrid this variable by
 *            default, if =-regrid_vars= was not set. Otherwise a
 *            variable is only regridded if both =-regrid_file= and
 *            =-regrid_vars= are set *and* the name of the variable is
 *            found in the set of names given with =-regrid_vars=.
 */
void Component::regrid(const std::string &module_name, IceModelVec *variable,
                       RegriddingFlag flag) {

  assert(variable != NULL);

  options::String regrid_file("-regrid_file", "regridding file name");

  options::StringSet regrid_vars("-regrid_vars",
                                 "comma-separated list of regridding variables",
                                 "");

  if (not regrid_file.is_set()) {
    return;
  }

  NCSpatialVariable &m = variable->metadata();

  if ((regrid_vars.is_set() and set_contains(regrid_vars, m.get_string("short_name"))) or
      (not regrid_vars.is_set() and flag == REGRID_WITHOUT_REGRID_VARS)) {

    verbPrintf(2, m_grid.com,
               "  %s: regridding '%s' from file '%s' ...\n",
               module_name.c_str(),
               m.get_string("short_name").c_str(), regrid_file->c_str());

    variable->regrid(regrid_file, CRITICAL);
  }
}

} // end of namespace pism
