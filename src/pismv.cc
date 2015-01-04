// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
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
"Ice sheet driver for PISM (SIA and SSA) verification.  Uses exact solutions\n"
"  to various coupled subsystems.  Computes difference between exact solution\n"
"  and numerical solution.  Can also just compute exact solution (-eo).\n"
"  Currently implements tests A, B, C, D, E, F, G, H, K, L.\n\n";

#include <cctype>               // toupper
#include <string>
#include <algorithm>            // std::transform()
#include <petscdmda.h>
#include "IceGrid.hh"
#include "verif/iceCompModel.hh"

#include "POConstant.hh"
#include "pism_options.hh"

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

// a wrapper that seems to be necessary to make std::transform below work
static inline char pism_toupper(char c)
{
    return (char)std::toupper(c);
}

int main(int argc, char *argv[]) {
  MPI_Comm com = MPI_COMM_WORLD;

  PetscInitializer petsc(argc, argv, help);
  com = PETSC_COMM_WORLD;
      
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();

    verbPrintf(2, com, "PISMV %s (verification mode)\n",
               PISM_Revision);

    if (options::Bool("-version", "stop after printing print PISM version")) {
      return 0;
    }

    std::string usage =
      "  pismv -test x [-no_report] [-eo] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -test x     SIA-type verification test (x = A|B|C|D|E|F|G|H|K|L)\n"
      "  -no_report  do not give error report at end of run\n"
      "  -eo         do not do numerical run; exact solution only\n"
      "(see User's Manual for tests I and J).\n";

    std::vector<std::string> required;
    required.push_back("-test");

    bool done = show_usage_check_req_opts(com, "pismv", required, usage);
    if (done) {
      return 0;
    }

    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    init_config(com, config, overrides, true);

    config.set_flag("use_eta_transformation", false);
    config.set_string("calendar", "none");

    IceGrid g(com, config);

    // determine test (and whether to report error)
    std::string testname = options::String("-test", "Specifies PISM verification test", "A");

    // transform to uppercase:
    transform(testname.begin(), testname.end(), testname.begin(), pism_toupper);

    // actually construct and run one of the derived classes of IceModel
    // run derived class for compensatory source SIA solutions
    // (i.e. compensatory accumulation or compensatory heating)
    IceCompModel m(g, config, overrides, testname[0]);
    m.setExecName("pismv");

    m.init();

    m.run();
    verbPrintf(2,com, "done with run\n");

    m.reportErrors();

    // provide a default output file name if no -o option is given.
    m.writeFiles("unnamed.nc");
  }
  catch (...) {
    handle_fatal_errors(com);
  }


  return 0;
}

