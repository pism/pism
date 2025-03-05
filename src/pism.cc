// Copyright (C) 2004--2024 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "util/io/IO_Flags.hh"
static char help[] =
  "Ice sheet driver for PISM ice sheet simulations, initialized from data.\n"
  "The basic PISM executable for evolution runs.\n";

#include <memory>
#include <petscsys.h>           // PETSC_COMM_WORLD

#include "pism/icemodel/IceModel.hh"
#include "pism/icemodel/IceEISModel.hh"
#include "pism/verification/iceCompModel.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"

#include "pism/util/Context.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/EnthalpyConverter.hh"

#include "pism/regional/IceRegionalModel.hh"

using namespace pism;

namespace eismint2 {

static void set_config_defaults(Config &config) {
  config.set_number("grid.Lx", 750); // in km
  config.set_number("grid.Ly", 750); // in km
  config.set_string("grid.periodicity", "none");
  config.set_string("grid.registration", "corner");
  config.set_string("stress_balance.sia.flow_law", "pb");

  config.set_string("energy.model", "none");

  // Set sea level elevation to -1e4 meters to remove ocean interaction
  config.set_number("sea_level.constant.value", -1e4);

  // purely SIA, and E=1
  config.set_number("stress_balance.sia.enhancement_factor", 1.0);

  // none use bed smoothing & bed roughness parameterization
  config.set_number("stress_balance.sia.bed_smoother.range", 0.0);

  // basal melt does not change computation of mass continuity or vertical velocity:
  config.set_flag("geometry.update.use_basal_melt_rate", false);

  // Make bedrock thermal material properties into ice properties.  Note that
  // zero thickness bedrock layer is the default, but we want the ice/rock
  // interface segment to have geothermal flux applied directly to ice without
  // jump in material properties at base.
  config.set_number("energy.bedrock_thermal.density",
                    config.get_number("constants.ice.density"));
  config.set_number("energy.bedrock_thermal.conductivity",
                    config.get_number("constants.ice.thermal_conductivity"));
  config.set_number("energy.bedrock_thermal.specific_heat_capacity",
                    config.get_number("constants.ice.specific_heat_capacity"));

  // no sliding + SIA
  config.set_string("stress_balance.model", "sia");
}
} // namespace eismint2

namespace verification {

//! Allocate the verification mode context. Uses ColdEnthalpyConverter.
std::shared_ptr<Context> context(MPI_Comm com, const std::string &prefix) {
  // unit system
  auto sys = std::make_shared<units::System>();

  // logger
  auto logger = logger_from_options(com);

  // configuration parameters
  auto config = config_from_options(com, *logger, sys);

  config->set_string("grid.periodicity", "none");
  config->set_string("grid.registration", "corner");

  set_config_from_options(sys, *config);
  config->resolve_filenames();

  print_config(*logger, 3, *config);

  auto time = std::make_shared<Time>(com, config, *logger, sys);

  auto EC = std::make_shared<ColdEnthalpyConverter>(*config);

  return std::make_shared<Context>(com, sys, config, EC, time, logger, prefix);
}

grid::Parameters grid_defaults(Config::Ptr config, char testname) {
  // This sets the defaults for each test; command-line options can override this.

  int Mx = 61, My = 61;
  double Lx = 1e3, Ly = 1e3;
  grid::Parameters P(*config, Mx, My, Lx, Ly);

  // use the cell corner grid registration
  P.registration = pism::grid::CELL_CORNER;
  // use the non-periodic grid:
  P.periodicity = pism::grid::NOT_PERIODIC;

  // equal spacing is the default for all the tests except K
  auto spacing    = pism::grid::EQUAL;
  double Lz       = config->get_number("grid.Lz");
  unsigned int Mz = config->get_number("grid.Mz");

  switch (testname) {
  case 'A':
  case 'B':
  case 'H':
    // use 2400km by 2400km by 4000m rectangular domain
    P.Lx = 1200e3;
    P.Ly = P.Lx;
    Lz   = 4000;
    break;
  case 'C':
  case 'D':
    // use 2000km by 2000km by 4000m rectangular domain
    P.Lx = 1000e3;
    P.Ly = P.Lx;
    Lz   = 4000;
    break;
  case 'F':
  case 'G':
  case 'L':
    // use 1800km by 1800km by 4000m rectangular domain
    P.Lx = 900e3;
    P.Ly = P.Lx;
    Lz   = 4000;
    break;
  case 'K':
  case 'O':
    // use 2000km by 2000km by 4000m rectangular domain, but make truely periodic
    config->set_number("grid.Mbz", 2);
    config->set_number("grid.Lbz", 1000);
    P.Lx          = 1000e3;
    P.Ly          = P.Lx;
    Lz            = 4000;
    P.periodicity = pism::grid::XY_PERIODIC;
    spacing       = pism::grid::QUADRATIC;
    break;
  case 'V':
    P.My          = 3;     // it's a flow-line setup
    P.Lx          = 500e3; // 500 km long
    P.periodicity = pism::grid::Y_PERIODIC;
    break;
  default:
    throw RuntimeError(PISM_ERROR_LOCATION, "desired test not implemented\n");
  }

  P.z = grid::compute_vertical_levels(Lz, Mz, spacing, config->get_number("grid.lambda"));
  return P;
}

std::shared_ptr<Grid> grid(std::shared_ptr<Context> ctx, char testname) {
  auto config = ctx->config();

  auto input_file_name = config->get_string("input.file");

  if (config->get_flag("input.bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION, "PISM does not support bootstrapping in verification mode");
  }

  if (not input_file_name.empty()) {
    auto r = grid::string_to_registration(config->get_string("grid.registration"));


    File input_file(ctx->com(), input_file_name, pism::io::PISM_NETCDF3, pism::io::PISM_READONLY);
    // get grid from a PISM input file
    return Grid::FromFile(ctx, input_file, { "enthalpy", "temp" }, r);
  }

  // use defaults set by grid_defaults()
  auto P = grid_defaults(config, testname);

  P.horizontal_size_and_extent_from_options(*config);
  P.vertical_grid_from_options(*config);
  P.ownership_ranges_from_options(*config, ctx->size());

  return std::make_shared<Grid>(ctx, P);
}
} // namespace verification

int main(int argc, char *argv[]) {

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  int exit_code = 0;
  try {
    // Note: EISMINT II experiments G and H are not supported.
    auto eisII = options::Keyword("-eisII",
                                  "EISMINT II experiment name",
                                  "A,B,C,D,E,F,I,J,K,L", "A");

    auto verification_test =
        options::Keyword("-test", "Specifies PISM verification test", "A,B,C,D,F,G,H,K,L,V", "A");

    if (eisII.is_set() and verification_test.is_set()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "please set -test or -eisII (not both)");
    }

    std::shared_ptr<Context> ctx;
    if (verification_test.is_set()) {
      ctx = verification::context(com, "pism");
    } else {
      ctx = context_from_options(com, "pism", false);
    }

    Logger::Ptr log = ctx->log();
    Config::Ptr config = ctx->config();

    std::vector<std::string> required_options{};
    if (eisII.is_set()) {
      // set defaults:
      eismint2::set_config_defaults(*config);

      // process -config_override
      DefaultConfig::Ptr overrides(new DefaultConfig(com,
                                                     "pism_overrides",
                                                     "-config_override",
                                                     ctx->unit_system()));
      overrides->init(*ctx->log());
      config->import_from(*overrides);
      // process command-line options
      set_config_from_options(ctx->unit_system(), *config);
      config->resolve_filenames();
    } else if (not verification_test.is_set()) {
      required_options.emplace_back("-i");
    }

    print_config(*ctx->log(), 3, *config);

    std::string usage =
      "  pism -i IN.nc [-bootstrap] [-regional] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i                         IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -bootstrap                 enable heuristics to produce an initial state from an incomplete input\n"
      "  -regional                  enable \"regional mode\"\n"
      "  -eisII [experiment]        enable EISMINT II mode\n"
      "  -test  [verification_test] enable verification mode\n"
      "notes:\n"
      "  * option -i is required\n"
      "  * if -bootstrap is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    {
      bool done = show_usage_check_req_opts(*log, "PISM (basic evolution run mode)" ,
                                            required_options, usage);
      if (done) {
        return 0;
      }
    }

    options::String profiling_log = options::String("-profile",
                                                    "Save detailed profiling data to a file.");

    if (profiling_log.is_set()) {
      ctx->profiling().start();
    }

    std::shared_ptr<Grid> grid;
    std::shared_ptr<IceModel> model;
    std::shared_ptr<IceCompModel> verification_model;

    if (verification_test.is_set()) {
      char test = verification_test.value()[0];
      grid = verification::grid(ctx, test);

      verification_model = std::make_shared<IceCompModel>(grid, ctx, test);
      model = verification_model;
    } else {
      grid = Grid::FromOptions(ctx);

      if (options::Bool("-regional", "enable regional (outlet glacier) mode")) {
        model = std::make_shared<IceRegionalModel>(grid, ctx);
      } else if (eisII.is_set()) {
        char experiment = eisII.value()[0];

        model = std::make_shared<IceEISModel>(grid, ctx, experiment);
      } else {
        model = std::make_shared<IceModel>(grid, ctx);
      }
    }

    model->init();

    auto list_type = options::Keyword("-list_diagnostics",
                                      "List available diagnostic quantities and stop.",
                                      "all,spatial,scalar,json",
                                      "all");

    if (list_type.is_set()) {
      model->list_diagnostics(list_type);
    } else {
      auto termination_reason = model->run();

      switch (termination_reason) {
      case PISM_CHEKPOINT:
        {
          exit_code = static_cast<int>(config->get_number("output.checkpoint.exit_code"));
          log->message(2, "... stopping (exit_code=%d) after saving the checkpoint file\n",
                       exit_code);
          break;
        }
      case PISM_SIGNAL:
        {
          exit_code = 0;
          break;
        }
      case PISM_DONE:
        {
          log->message(2, "... done with the run\n");
          model->save_results();
          exit_code = 0;

          if (verification_model and
              not options::Bool("-no_report", "do not print the error report")) {
            verification_model->reportErrors();
          }
          break;
        }
      }
    }
    print_unused_parameters(*log, 3, *config);

    if (profiling_log.is_set()) {
      ctx->profiling().report(profiling_log);
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
