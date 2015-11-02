// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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

//! @file pismo.cc A regional (outlet glacier) model form of PISM.
/*! @file pismo.cc
The classes used by this driver program modify basic PISM whole ice sheet modeling assumptions.
Normally in PISM the ice sheet occupies a continent which is surrounded by
ocean.  Or at least PISM assumes that the edge of the computational domain is in
a region with strong ablation that the ice will not cross.

Here, by contrast, we add a strip around the edge of the computational domain
(variable `no_model_mask` and option `-no_model_strip`).  Various
simplifications and boundary conditions are enforced in this script:
* the surface gradient computation is made trivial,
* the driving stress does not change during the run but instead comes from
the gradient of a saved surface elevation, and
* the base is made strong so that no sliding occurs.

Also options `-force_to_thickness_file` and variable `ftt_mask` play a role in isolating
the modeled outlet glacier.  But there is no code here for that purpose.
Instead see the ForceThickness surface model modifier class.
 */

static char help[] =
  "Ice sheet driver for PISM regional (outlet glacier) simulations, initialized\n"
  "from data.\n";

#include <petscsys.h>

#include "base/util/IceGrid.hh"
#include "IceRegionalModel.hh"

#include "base/util/PISMConfig.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_options.hh"
#include "base/util/Context.hh"

int main(int argc, char *argv[]) {

  using namespace pism;

  MPI_Comm com = MPI_COMM_WORLD;

  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();

    verbPrintf(2,com, "PISMO %s (regional outlet-glacier run mode)\n",
               PISM_Revision);

    if (options::Bool("-version", "stop after printing print PISM version")) {
      return 0;
    }

    bool input_file_set = options::Bool("-i", "input file name");
    std::string usage =
      "  pismo -i IN.nc [-bootstrap] [-no_model_strip X] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i IN.nc   is input file in NetCDF format: contains PISM-written model state\n"
      "  -bootstrap enables heuristics used to produce an initial state from incomplete input\n"
      "  -no_model_strip X (re-)set width of no-model strip along edge of\n"
      "              computational domain to X km\n"
      "notes:\n"
      "  * option -i is required\n"
      "  * if -bootstrap is used then also '-Mx A -My B -Mz C -Lz D' are required\n";

    if (not input_file_set) {
      throw RuntimeError("options -i is required\n\n" + usage);
    }

    bool done = show_usage_check_req_opts(com, "pismo",
                                          std::vector<std::string>(), // no required options
                                          usage);
    if (done) {
      return 0;
    }

    Context::Ptr ctx = context_from_options(com, "pismo");
    Config::Ptr config = ctx->config();

    // initialize the ice dynamics model
    IceGrid::Ptr g = IceGrid::FromOptions(ctx);
    IceRegionalModel m(g, ctx);

    m.init();

    m.run();

    verbPrintf(2,com, "... done with run\n");

    // provide a default output file name if no -o option is given.
    m.writeFiles("unnamed_regional.nc");

    print_unused_parameters(*ctx->log(), 3, *config);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
