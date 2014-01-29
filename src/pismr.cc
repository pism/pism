// Copyright (C) 2004-2011, 2013, 2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

static char help[] =
  "Ice sheet driver for PISM ice sheet simulations, initialized from data.\n"
  "The basic PISM executable for evolution runs.\n";

#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"

#include "pism_options.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMR %s (basic evolution run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
    std::string usage =
      "  pismr {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismr", usage); CHKERRQ(ierr);
    } else {
      std::vector<std::string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pismr", required, usage); CHKERRQ(ierr);
    }

    PISMUnitSystem unit_system(NULL);
    PISMConfig config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    ierr = init_config(com, config, overrides, true); CHKERRQ(ierr);

    IceGrid g(com, config);
    IceModel m(g, config, overrides);

    ierr = m.setExecName("pismr"); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);
    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamed.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
