// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016, 2017 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <mpi.h>
extern "C"{
#include "cdipio.h"
#include "yaxt.h"
#include "cdi.h"
}
#endif

#include "pism/util/cdipio/CDIPIOInitializer.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);

  MPI_Comm com = MPI_COMM_WORLD;
  MPI_Comm local_comm;
  int nwriters, IOmode;
  bool async;
  {
  petsc::Initializer petsc(argc, argv, help);
  {
      Context::Ptr ctx = initial_context_from_options(com, "pismr");
      nwriters = ctx->get_n_writers();
      IOmode = ctx->get_IOmode();
      async = ctx->get_async();
  }
  }
  {
  cdipio::Initializer cdipio(nwriters, IOmode, com, async);
  local_comm = cdipio.get_comp_comm();

  try {
#if (Pism_USE_CDIPIO==1)
    if (local_comm != MPI_COMM_NULL) {
      cdipio.activate_namespace();
      PETSC_COMM_WORLD = local_comm;
      petsc::Initializer petsc(argc, argv, help);
#else
    local_comm = PETSC_COMM_WORLD;
#endif
    Context::Ptr ctx = context_from_options(local_comm, "pismr");
    Logger::Ptr log = ctx->log();

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

    options::String profiling_log = options::String("-profile",
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
#if (Pism_USE_CDIPIO==1)
   }
#endif
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }
  } // Finalize YAXT and CDI-PIO
  MPI_Finalize();
  return 0;
}
