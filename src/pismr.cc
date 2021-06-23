// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016, 2017, 2021 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism/util/IceGrid.hh"
#include "pism/icemodel/IceModel.hh"
#include "pism/util/Config.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Context.hh"
#include "pism/util/Profiling.hh"

#include "pism/regional/IceGrid_Regional.hh"
#include "pism/regional/IceRegionalModel.hh"

#if (Pism_USE_CDIPIO==1)
#include "pism/util/cdipio/CDIPIOInitializer.hh"
#endif

using namespace pism;

int main(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);

  // The communicator for processes involved in computation (as opposed to I/O).
  MPI_Comm com = MPI_COMM_WORLD;
  {
#if (Pism_USE_CDIPIO==1)
    MPI_Comm world = MPI_COMM_WORLD;
    int cdipio_nwriters;
    std::string cdipio_io_mode;
    {
      // Initialize PETSc on MPI_COMM_WORLD to be able to use PetscOptionsXXX in the code
      // below. Will be finalized at the end of this code block.
      petsc::Initializer petsc(argc, argv, help);

      auto log = logger_from_options(com);
      units::System::Ptr sys(new units::System);
      auto config = config_from_options(world, *log, sys);

      cdipio_nwriters = config->get_number("output.cdi_pio.n_writers");
      cdipio_io_mode  = config->get_string("output.cdi_pio.mode");
    }
    std::unique_ptr<cdipio::Initializer> cdipio;
    cdipio.reset(new cdipio::Initializer(cdipio_nwriters, cdipio_io_mode, world));
    com = cdipio->comp_comm();
    if (com == MPI_COMM_NULL) {
      // com is null if this process is a part of the I/O sub-communicator
      // YAXT and CDI-PIO are finalized here
      cdipio.reset(nullptr);
      MPI_Finalize();
      return 0;
    }
    cdipio->activate_namespace();
#endif
    PETSC_COMM_WORLD = com;
    std::unique_ptr<petsc::Initializer> petsc;
    petsc.reset(new petsc::Initializer(argc, argv, help));
    try {
      auto ctx = context_from_options(com, "pismr");
      auto log = ctx->log();

      std::string usage =
        "  pismr -i IN.nc [-bootstrap] [-regional] [OTHER PISM & PETSc OPTIONS]\n"
        "where:\n"
        "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
        "  -bootstrap  enable heuristics to produce an initial state from an incomplete input\n"
        "  -regional   enable \"regional mode\"\n"
        "notes:\n"
        "  * option -i is required\n"
        "  * if -bootstrap is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
      {
        std::vector<std::string> required(1, "-i");

        bool done = show_usage_check_req_opts(*log, "PISMR (basic evolution run mode)" ,
                                              required, usage);
        if (done) {
          return 0;
        }
      }

      auto profiling_log = options::String("-profile",
                                           "Save detailed profiling data to a file.");

      Config::Ptr config = ctx->config();

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
        model.reset(new IceModel(grid, ctx));
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
        model->run();

        log->message(2, "... done with run\n");

        model->save_results();
        model->close_files();
      }
      print_unused_parameters(*log, 3, *config);

      if (profiling_log.is_set()) {
        ctx->profiling().report(profiling_log);
      }
    }
    catch (...) {
      handle_fatal_errors(com);
      return 1;
    }
    // Enforce the order: PETSc gets finalized first, then CDI-PIO and YAXT, then MPI.
    petsc.reset(nullptr);
#if (Pism_USE_CDIPIO==1)
    cdipio.reset(nullptr);
#endif
  }
  MPI_Finalize();
  return 0;
}
