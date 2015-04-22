// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "Ice sheet driver for EISMINT II, and other constant climate, simplified geometry\n"
  "intercomparison simulations.\n";

#include <petscsys.h>

#include "base/iceModel.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfig.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_options.hh"
#include "eismint/iceEISModel.hh"

using namespace pism;

int main(int argc, char *argv[]) {

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    verbosityLevelFromOptions();

    verbPrintf(2,com, "PISMS %s (simplified geometry mode)\n",
               PISM_Revision);

    if (options::Bool("-version", "stop after printing print PISM version")) {
      return 0;
    }

    std::string usage =
      "  pisms [-eisII x] [OTHER PISM & PETSc OPTIONS]\n"
      "where major option chooses type of simplified experiment:\n"
      "  -eisII x    choose EISMINT II experiment (x = A|B|C|D|E|F|G|H|I|J|K|L)\n";

    std::vector<std::string> required;
    required.clear(); // no actually required options; "-eisII A" is default

    bool done = show_usage_check_req_opts(com, "pisms", required, usage);
    if (done) {
      return 0;
    }

    units::System::Ptr unit_system(new units::System);
    DefaultConfig::Ptr config(new DefaultConfig(com, "pism_config", "-config", unit_system)),
      overrides(new DefaultConfig(com, "pism_overrides", "-config_override", unit_system));
    overrides->init();
    config->init_with_default();
    config->import_from(*overrides);
    config->set_from_options();
    print_config(3, com, *config);

    config->set_string("calendar", "none");

    IceGrid g(com, config);
    IceEISModel m(g, config, overrides);

    m.init();

    m.run();

    verbPrintf(2,com, "... done with run \n");

    // provide a default output file name if no -o option is given.
    m.writeFiles("unnamed.nc");

    print_unused_parameters(3, com, *config);
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
