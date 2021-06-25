// Copyright (C) 2004-2017, 2019, 2020 Jed Brown, Ed Bueler and Constantine Khroulev
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
"  and numerical solution.\n"
"  Currently implements tests A, B, C, D, E, F, G, H, K, L.\n\n";

#include <string>

#include "pism/util/IceGrid.hh"
#include "pism/util/Config.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/verification/iceCompModel.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Time.hh"
#include "pism/util/EnthalpyConverter.hh"

using namespace pism;

//! Allocate the PISMV (verification) context. Uses ColdEnthalpyConverter.
std::shared_ptr<Context> pismv_context(MPI_Comm com, const std::string &prefix) {
  // unit system
  units::System::Ptr sys(new units::System);

  // logger
  Logger::Ptr logger = logger_from_options(com);

  // configuration parameters
  Config::Ptr config = config_from_options(com, *logger, sys);

  config->set_string("grid.periodicity", "none");
  config->set_string("grid.registration", "corner");

  set_config_from_options(sys, *config);

  print_config(*logger, 3, *config);

  Time::Ptr time = time_from_options(com, config, sys);

  EnthalpyConverter::Ptr EC = EnthalpyConverter::Ptr(new ColdEnthalpyConverter(*config));

  return std::shared_ptr<Context>(new Context(com, sys, config, EC, time, logger, prefix));
}

GridParameters pismv_grid_defaults(Config::Ptr config,
                                   char testname) {
  // This sets the defaults for each test; command-line options can override this.

  GridParameters P;

  // use the cell corner grid registration
  P.registration = CELL_CORNER;
  // use the non-periodic grid:
  P.periodicity = NOT_PERIODIC;
  // equal spacing is the default for all the tests except K
  P.Lx = config->get_number("grid.Lx");
  P.Ly = config->get_number("grid.Ly");

  P.Mx = config->get_number("grid.Mx");
  P.My = config->get_number("grid.My");

  SpacingType spacing = EQUAL;
  double Lz = config->get_number("grid.Lz");
  unsigned int Mz = config->get_number("grid.Mz");

  switch (testname) {
  case 'A':
  case 'B':
  case 'H':
    // use 2400km by 2400km by 4000m rectangular domain
    P.Lx = 1200e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'C':
  case 'D':
    // use 2000km by 2000km by 4000m rectangular domain
    P.Lx = 1000e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'F':
  case 'G':
  case 'L':
    // use 1800km by 1800km by 4000m rectangular domain
    P.Lx = 900e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'K':
  case 'O':
    // use 2000km by 2000km by 4000m rectangular domain, but make truely periodic
    config->set_number("grid.Mbz", 2);
    config->set_number("grid.Lbz", 1000);
    P.Lx = 1000e3;
    P.Ly = P.Lx;
    Lz = 4000;
    P.periodicity = XY_PERIODIC;
    spacing = QUADRATIC;
    break;
  case 'V':
    P.My = 3;             // it's a flow-line setup
    P.Lx = 500e3;            // 500 km long
    P.periodicity = Y_PERIODIC;
    break;
  default:
    throw RuntimeError(PISM_ERROR_LOCATION, "desired test not implemented\n");
  }

  P.z = IceGrid::compute_vertical_levels(Lz, Mz, spacing,
                                         config->get_number("grid.lambda"));
  return P;
}

IceGrid::Ptr pismv_grid(std::shared_ptr<Context> ctx, char testname) {
  auto config = ctx->config();

  auto input_file = config->get_string("input.file");

  if (config->get_flag("input.bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION, "pismv does not support bootstrapping");
  }

  if (not input_file.empty()) {
    GridRegistration r = string_to_registration(ctx->config()->get_string("grid.registration"));

    // get grid from a PISM input file
    return IceGrid::FromFile(ctx, input_file, {"enthalpy", "temp"}, r);
  } else {
    // use defaults set by pismv_grid_defaults()
    GridParameters P = pismv_grid_defaults(ctx->config(), testname);
    P.horizontal_size_from_options();
    P.horizontal_extent_from_options(ctx->unit_system());
    P.vertical_grid_from_options(ctx->config());
    P.ownership_ranges_from_options(ctx->size());

    return IceGrid::Ptr(new IceGrid(ctx, P));
  }
}

int main(int argc, char *argv[]) {
  MPI_Comm com = MPI_COMM_WORLD;

  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;
      
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    std::shared_ptr<Context> ctx = pismv_context(com, "pismv");
    Logger::Ptr log = ctx->log();

    std::string usage =
      "  pismv -test x [-no_report] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -test x     SIA-type verification test (x = A|B|C|D|F|G|H|K|L)\n"
      "  -no_report  do not give error report at end of run\n"
      "(see User's Manual for tests I and J).\n";

    std::vector<std::string> required(1, "-test");

    bool done = show_usage_check_req_opts(*log, "PISMV (verification mode)", required, usage);
    if (done) {
      return 0;
    }

    Config::Ptr config = ctx->config();

    // determine test (and whether to report error)
    std::string testname = options::Keyword("-test", "Specifies PISM verification test",
                                            "A,B,C,D,F,G,H,K,L,V", "A");

    IceGrid::Ptr g = pismv_grid(ctx, testname[0]);

    IceCompModel m(g, ctx, testname[0]);

    m.init();

    m.run();
    log->message(2, "done with run\n");

    if (not options::Bool("-no_report", "do not print the error report")) {
      m.reportErrors();
    }

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

