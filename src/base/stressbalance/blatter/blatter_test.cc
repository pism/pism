// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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
  "The executable for testing the Blatter stress balance solver.\n";

#include "base/stressbalance/blatter/BlatterStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/io/PIO.hh"
#include "base/util/PISMVars.hh"
#include "base/rheology/FlowLaw.hh"
#include "base/enthalpyConverter.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMTime.hh"
#include "coupler/ocean/POConstant.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/io_helpers.hh"


#include "base/util/petscwrappers/PetscInitializer.hh"

using namespace pism;

//! \brief Read data from an input file.
static void read_input_data(const std::string &filename,
                            Vars &variables,
                            EnthalpyConverter::Ptr EC) {
  // Get the names of all the variables allocated earlier:
  std::set<std::string> vars = variables.keys();

  for (std::set<std::string>::iterator i = vars.begin(); i != vars.end(); ++i) {
    if (*i != "enthalpy") {
      variables.get_shared(*i)->regrid(filename, CRITICAL);
    } else {
      const double
        T = 263.15,
        omega = 0.0,
        pressure = 0.0;
      variables.get_shared(*i)->set(EC->enthalpy(T, omega, pressure));
    }
  }
}

//! \brief Write data to an output file.
static void write_input_fields(std::string filename, Vars &variables) {
  // Get the names of all the variables allocated earlier:
  std::set<std::string> vars = variables.keys();

  for (std::set<std::string>::iterator i = vars.begin(); i != vars.end(); ++i) {
    variables.get(*i)->write(filename);
  }
}

//! \brief Allocate IceModelVec2S variables.
static void allocate_variables(IceGrid::ConstPtr grid, Vars &variables) {

  IceModelVec2S::Ptr thk(new IceModelVec2S);
  IceModelVec2S::Ptr topg(new IceModelVec2S);
  IceModelVec2S::Ptr tauc(new IceModelVec2S);
  IceModelVec3::Ptr enthalpy(new IceModelVec3);

  thk->create(grid, "thk", WITH_GHOSTS);
  thk->set_attrs("", "ice thickness",
                 "m", "land_ice_thickness");
  variables.add_shared(thk);

  topg->create(grid, "topg", WITH_GHOSTS);
  topg->set_attrs("", "bedrock surface elevation",
                  "m", "bedrock_altitude");
  variables.add_shared(topg);

  tauc->create(grid, "tauc", WITH_GHOSTS);
  tauc->set_attrs("diagnostic",
                  "yield stress for basal till (plastic or pseudo-plastic model)",
                  "Pa", "");
  variables.add_shared(tauc);

  enthalpy->create(grid, "enthalpy", WITH_GHOSTS, 1);
  enthalpy->set_attrs("model_state",
                      "ice enthalpy (includes sensible heat, latent heat, pressure)",
                      "J kg-1", "");
  variables.add_shared(enthalpy);
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank,size
  PetscLogStage cold, hot;

  petsc::Initializer petsc(argc, argv, help);

  try {
    ierr = PetscLogStageRegister("Cold", &cold); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("Hot ", &hot); CHKERRQ(ierr);

    com = PETSC_COMM_WORLD;

    /* This explicit scoping forces destructors to be called before PetscFinalize() */
    Context::Ptr ctx = context_from_options(com, "blatter_test");
    Config::Ptr config = ctx->config();
    EnthalpyConverter::Ptr EC = ctx->enthalpy_converter();


    ctx->log()->message(2,
                        "BLATTER_TEST: testing the Blatter stress balance solver.\n");

    if (options::Bool("-usage", "print usage") or
        options::Bool("-help", "print usage")) {
      ctx->log()->message(2,
                          "\n"
                          "usage of BLATTER_TEST:\n"
                          "  run blatter_test -i input.nc -o output.nc\n"
                          "\n");
    }

    options::String input_file("-i", "Set the input file name");
    if (not input_file.is_set()) {
      throw RuntimeError("BLATTER_TEST ERROR: -i is required.");
    }

    options::String output_file("-o", "Set the output file name", "blatter_test.nc");

    IceGrid::Ptr grid = IceGrid::FromOptions(ctx);

    allocate_variables(grid, grid->variables());

    ocean::Constant ocean(grid);

    read_input_data(input_file, grid->variables(), EC);

    IceModelVec2S melange_back_pressure;
    melange_back_pressure.create(grid, "melange_back_pressure", WITHOUT_GHOSTS);
    melange_back_pressure.set_attrs("boundary_condition",
                                    "melange back pressure fraction",
                                    "1", "");
    melange_back_pressure.set(0.0);

    PetscLogStagePush(cold);
    stressbalance::BlatterStressBalance blatter(grid, EC);
    // Initialize the Blatter solver:
    blatter.init();
    blatter.update(false, melange_back_pressure);
    PetscLogStagePop();

    if (options::Bool("-compare", "Compare 'cold' and 'hot' runs.")) {
      PetscLogStagePush(hot);
      blatter.update(false, melange_back_pressure);
      PetscLogStagePop();
    }

    // Write results to an output file:
    PIO pio(grid->com, "netcdf3");

    pio.open(output_file, PISM_READWRITE_MOVE);
    io::define_time(pio, config->get_string("time_dimension_name"),
                    grid->ctx()->time()->calendar(),
                    grid->ctx()->time()->CF_units_string(),
                    ctx->unit_system());
    io::append_time(pio, config->get_string("time_dimension_name"), 0.0);
    // Ensure that time was appended (io::append_time above schedules a write, but does not perform
    // any I/O).
    pio.close();

    // input fields (bed elevation, ice thickness, basal yield stress, enthalpy)
    write_input_fields(output_file, grid->variables());

    // 3D velocity components on the regular grid
    blatter.velocity_u().write(output_file);
    blatter.velocity_v().write(output_file);
    // 3D velocity components on the sigma grid
    blatter.velocity_u_sigma().write(output_file);
    blatter.velocity_v_sigma().write(output_file);
    // and the vertically-averaged velocity
    blatter.velocity().write(output_file);
  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
