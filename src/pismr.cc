// Copyright (C) 2004-2011, 2013, 2014, 2015 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  PetscErrorCode ierr;

  MPI_Comm com = MPI_COMM_WORLD;
  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    verbosityLevelFromOptions();

    verbPrintf(2,com, "PISMR %s (basic evolution run mode)\n",
               PISM_Revision);
    stop_on_version_option();

    bool iset = options::Bool("-i", "input file name");
    bool bfset = options::Bool("-boot_file", "bootstrapping file name");
    std::string usage =
      "  pismr {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((iset == false) && (bfset == false)) {
      ierr = PetscPrintf(com,
                         "\nPISM ERROR: one of options -i,-boot_file is required\n\n");
      PISM_PETSC_CHK(ierr, "PetscPrintf");
      show_usage_and_quit(com, "pismr", usage);
    } else {
      std::vector<std::string> required;  required.clear();
      show_usage_check_req_opts(com, "pismr", required, usage);
    }

    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    init_config(com, config, overrides, true);

    IceGrid g(com, config);
    IceModel m(g, config, overrides);

    m.setExecName("pismr");

    m.init();

    bool print_list_and_stop = options::Bool("-list_diagnostics",
                                             "List available diagnostic quantities and stop");

    if (print_list_and_stop) {
      m.list_diagnostics();
    } else {
      m.run();

      verbPrintf(2,com, "... done with run\n");
      // provide a default output file name if no -o option is given.
      m.writeFiles("unnamed.nc");
    }
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
