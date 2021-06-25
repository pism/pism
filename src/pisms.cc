// Copyright (C) 2004-2019, 2021 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism/icemodel/IceModel.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Config.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/icemodel/IceEISModel.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Time.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/io/File.hh"

using namespace pism;

static void pisms_config_defaults(Config &config) {
  config.set_number("grid.Lx", 750e3);
  config.set_number("grid.Ly", 750e3);
  config.set_string("grid.periodicity", "none");
  config.set_string("grid.registration", "corner");
  config.set_string("stress_balance.sia.flow_law", "pb");
}

int main(int argc, char *argv[]) {

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "");
    pisms_config_defaults(*ctx->config());
    Logger::Ptr log = ctx->log();

    std::string usage =
      "  pisms [-eisII x] [OTHER PISM & PETSc OPTIONS]\n"
      "where major option chooses type of simplified experiment:\n"
      "  -eisII x    choose EISMINT II experiment (x = A|B|C|D|E|F|I|J|K|L)\n";

    bool done = show_usage_check_req_opts(*log, "PISMS (simplified geometry mode)",
                                          std::vector<std::string>(), // no required options
                                          usage);
    if (done) {
      return 0;
    }

    // Note: experiments G and H are not supported.
    std::string experiment = options::Keyword("-eisII", "EISMINT II experiment name",
                                              "A,B,C,D,E,F,I,J,K,L", "A");

    Config::Ptr config = ctx->config();

    IceGrid::Ptr g = IceGrid::FromOptions(ctx);
    IceEISModel m(g, ctx, experiment[0]);

    m.init();

    m.run();

    log->message(2, "... done with run \n");

    // provide a default output file name if no -o option is given.
    m.save_results();

    print_unused_parameters(*log, 3, *config);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
