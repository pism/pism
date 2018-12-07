// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Ed Bueler and Constantine Khroulev
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
  "Tests BedThermalUnit using Test K, without IceModel.\n\n";

#include "pism/util/pism_options.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/verification/BTU_Verification.hh"
#include "pism/energy/BTU_Minimal.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/ConfigInterface.hh"

#include "pism/verification/tests/exactTestK.h"

#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Context.hh"
#include "pism/util/Config.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Logger.hh"

//! Allocate the PISMV (verification) context. Uses ColdEnthalpyConverter.
pism::Context::Ptr btutest_context(MPI_Comm com, const std::string &prefix) {
  using namespace pism;

  // unit system
  units::System::Ptr sys(new units::System);

  // logger
  Logger::Ptr logger = logger_from_options(com);

  // configuration parameters
  Config::Ptr config = config_from_options(com, *logger, sys);

  config->set_string("time.calendar", "none");
  // when IceGrid constructor is called, these settings are used
  config->set_double("time.start_year", 0.0);
  config->set_double("time.run_length", 1.0);

  set_config_from_options(*config);

  print_config(*logger, 3, *config);

  Time::Ptr time = time_from_options(com, config, sys);

  EnthalpyConverter::Ptr EC = EnthalpyConverter::Ptr(new ColdEnthalpyConverter(*config));

  return Context::Ptr(new Context(com, sys, config, EC, time, logger, prefix));
}

int main(int argc, char *argv[]) {

  using namespace pism;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    Context::Ptr ctx = btutest_context(com, "btutest");
    Logger::Ptr log = ctx->log();

    std::string usage =
      "  btutest -Mbz NN -Lbz 1000.0 [-o OUT.nc -ys A -ye B -dt C -Mz D -Lz E]\n"
      "where these are required because they are used in BedThermalUnit:\n"
      "  -Mbz           number of bedrock thermal layer levels to use\n"
      "  -Lbz 1000.0    depth of bedrock thermal layer (required; Lbz=1000.0 m in Test K)\n"
      "and these are allowed:\n"
      "  -o             output file name; NetCDF format\n"
      "  -ys            start year in using Test K\n"
      "  -ye            end year in using Test K\n"
      "  -dt            time step B (= positive float) in years\n"
      "  -Mz            number of ice levels to use\n"
      "  -Lz            height of ice/atmospher box\n";

    // check required options
    std::vector<std::string> required(1, "-Mbz");

    bool done = show_usage_check_req_opts(*log, "BTUTEST %s (test program for BedThermalUnit)",
                                          required, usage);
    if (done) {
      return 0;
    }

    Config::Ptr config = ctx->config();

    // Mbz and Lbz are used by the BedThermalUnit, not by IceGrid
    config->set_double("grid.Mbz", 11);
    config->set_double("grid.Lbz", 1000);

    // Set default vertical grid parameters.
    config->set_double("grid.Mz", 41);
    config->set_double("grid.Lz", 4000);

    log->message(2,
                 "  initializing IceGrid from options ...\n");

    GridParameters P(config);
    P.Mx = 3;
    P.My = P.Mx;
    P.Lx = 1500e3;
    P.Ly = P.Lx;

    P.vertical_grid_from_options(config);
    P.ownership_ranges_from_options(ctx->size());

    // create grid and set defaults
    IceGrid::Ptr grid(new IceGrid(ctx, P));

    ctx->time()->init(*log);

    auto outname = config->get_string("output.file_name");

    options::Real dt_years("-dt", "Time-step, in years", 1.0);

    // allocate tools and IceModelVecs
    IceModelVec2S bedtoptemp, heat_flux_at_ice_base;
    {
      heat_flux_at_ice_base.create(grid, "upward_heat_flux_at_ice_base", WITHOUT_GHOSTS);
      heat_flux_at_ice_base.set_attrs("",
                     "upward geothermal flux at bedrock thermal layer base",
                     "W m-2", "");
      heat_flux_at_ice_base.metadata().set_string("glaciological_units", "mW m-2");

      bedtoptemp.create(grid, "bedtoptemp", WITHOUT_GHOSTS);
      bedtoptemp.set_attrs("",
                            "temperature at top of bedrock thermal layer",
                            "K", "");
    }

    // initialize BTU object:
    energy::BTUGrid bedrock_grid = energy::BTUGrid::FromOptions(ctx);

    energy::BedThermalUnit::Ptr btu;

    if (bedrock_grid.Mbz > 1) {
      btu.reset(new energy::BTU_Verification(grid, bedrock_grid, 'K', false));
    } else {
      btu.reset(new energy::BTU_Minimal(grid));
    }

    InputOptions opts = process_input_options(com, config);
    btu->init(opts);

    double dt_seconds = units::convert(ctx->unit_system(), dt_years, "years", "seconds");

    // worry about time step
    int  N = (int)ceil((ctx->time()->end() - ctx->time()->start()) / dt_seconds);
    dt_seconds = (ctx->time()->end() - ctx->time()->start()) / (double)N;
    log->message(2,
                 "  user set timestep of %.4f years ...\n"
                 "  reset to %.4f years to get integer number of steps ... \n",
                 dt_years.value(), units::convert(ctx->unit_system(), dt_seconds, "seconds", "years"));
    MaxTimestep max_dt = btu->max_timestep(0.0);
    log->message(2,
                 "  BedThermalUnit reports max timestep of %.4f years ...\n",
                 units::convert(ctx->unit_system(), max_dt.value(), "seconds", "years"));

    // actually do the time-stepping
    log->message(2,"  running ...\n");
    for (int n = 0; n < N; n++) {
      // time at start of time-step
      const double time = ctx->time()->start() + dt_seconds * (double)n;

      // compute exact ice temperature at z=0 at time y
      {
        IceModelVec::AccessList list(bedtoptemp);
        for (Points p(*grid); p; p.next()) {
          const int i = p.i(), j = p.j();

          bedtoptemp(i,j) = exactK(time, 0.0, 0).T;
        }
      }
      // no need to update ghost values

      // update the temperature inside the thermal layer using bedtoptemp
      btu->update(bedtoptemp, time, dt_seconds);
      log->message(2,".");
    }

    log->message(2, "\n  done ...\n");

    // compute final output heat flux G_0 at z=0
    heat_flux_at_ice_base.copy_from(btu->flux_through_top_surface());

    // get, and tell stdout, the correct answer from Test K
    const double FF = exactK(ctx->time()->end(), 0.0, 0).F;
    log->message(2,
                 "  exact Test K reports upward heat flux at z=0, at end time %s, as G_0 = %.7f W m-2;\n",
                 ctx->time()->end_date().c_str(), FF);

    // compute numerical error
    heat_flux_at_ice_base.shift(-FF);
    double max_error = heat_flux_at_ice_base.norm(NORM_INFINITY);
    double avg_error = heat_flux_at_ice_base.norm(NORM_1);
    heat_flux_at_ice_base.shift(+FF); // shift it back for writing
    avg_error /= (grid->Mx() * grid->My());
    log->message(2,
                 "case dt = %.5f:\n", dt_years.value());
    log->message(1,
                 "NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:\n");
    log->message(1,
                 "bheatflx0  :       max    prcntmax          av\n");
    log->message(1,
                 "           %11.7f  %11.7f  %11.7f\n",
                 max_error, 100.0*max_error/FF, avg_error);
    log->message(1, "NUM ERRORS DONE\n");

    PIO file(grid->com, grid->ctx()->config()->get_string("output.format"),
            outname, PISM_READWRITE_MOVE);

    io::define_time(file, *ctx);
    io::append_time(file, *ctx->config(), ctx->time()->current());

    btu->write_model_state(file);

    bedtoptemp.write(file);
    heat_flux_at_ice_base.write(file);

    log->message(2, "done.\n");
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
