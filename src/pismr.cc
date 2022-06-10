// Copyright (C) 2004--2022 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <memory>
#include <petscsys.h>           // PETSC_COMM_WORLD

#include "pism/icemodel/IceModel.hh"
#include "pism/icemodel/IceEISModel.hh"
#include "pism/util/Config.hh"
#include "pism/util/IceGrid.hh"

#include "pism/util/Context.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"

#include "pism/regional/IceGrid_Regional.hh"
#include "pism/regional/IceRegionalModel.hh"

using namespace pism;

static void set_eismint2_config_defaults(Config &config) {
  config.set_number("grid.Lx", 750e3);
  config.set_number("grid.Ly", 750e3);
  config.set_string("grid.periodicity", "none");
  config.set_string("grid.registration", "corner");
  config.set_string("stress_balance.sia.flow_law", "pb");

  config.set_flag("energy.temperature_based", true);

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

    std::shared_ptr<Context> ctx = context_from_options(com, "pismr", false);

    Logger::Ptr log = ctx->log();
    Config::Ptr config = ctx->config();

    std::vector<std::string> required_options{};
    if (eisII.is_set()) {
      // set defaults:
      set_eismint2_config_defaults(*config);

      // process -config_override
      DefaultConfig::Ptr overrides(new DefaultConfig(com,
                                                     "pism_overrides",
                                                     "-config_override",
                                                     ctx->unit_system()));
      overrides->init(*ctx->log());
      config->import_from(*overrides);
      // process command-line options
      set_config_from_options(ctx->unit_system(), *config);
    } else {
      required_options.emplace_back("-i");
    }

    print_config(*ctx->log(), 3, *config);

    std::string usage =
      "  pismr -i IN.nc [-bootstrap] [-regional] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i                   IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -bootstrap           enable heuristics to produce an initial state from an incomplete input\n"
      "  -regional            enable \"regional mode\"\n"
      "  -eisII [experiment]  enable EISMINT II mode\n"
      "notes:\n"
      "  * option -i is required\n"
      "  * if -bootstrap is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    {
      bool done = show_usage_check_req_opts(*log, "PISMR (basic evolution run mode)" ,
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

    IceGrid::Ptr grid;
    std::unique_ptr<IceModel> model;

    if (options::Bool("-regional", "enable regional (outlet glacier) mode")) {
      grid = regional_grid_from_options(ctx);
      model.reset(new IceRegionalModel(grid, ctx));
    } else {
      grid = IceGrid::FromOptions(ctx);

      if (eisII.is_set()) {
        char experiment = eisII.value()[0];
        model.reset(new IceEISModel(grid, ctx, experiment));
      } else {
        model.reset(new IceModel(grid, ctx));
      }
    }

    model->init();

    const bool
      list_ascii = options::Bool("-list_diagnostics",
                                 "List available diagnostic quantities and stop"),
      list_json = options::Bool("-list_diagnostics_json",
                                "List available diagnostic quantities (JSON format) and stop");

    if (list_ascii) {
      model->list_diagnostics();
    } else if (list_json) {
      model->list_diagnostics_json();
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
