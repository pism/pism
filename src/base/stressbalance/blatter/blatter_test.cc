// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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
#include "base/rheology/flowlaws.hh"
#include "base/enthalpyConverter.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMTime.hh"
#include "coupler/ocean/POConstant.hh"
#include "base/util/error_handling.hh"

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
static void write_data(std::string filename, Vars &variables) {
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
  variables.add(*thk);

  topg->create(grid, "topg", WITH_GHOSTS);
  topg->set_attrs("", "bedrock surface elevation",
                  "m", "bedrock_altitude");
  variables.add(*topg);

  tauc->create(grid, "tauc", WITH_GHOSTS);
  tauc->set_attrs("diagnostic",
                  "yield stress for basal till (plastic or pseudo-plastic model)",
                  "Pa", "");
  variables.add(*tauc);

  enthalpy->create(grid, "enthalpy", WITH_GHOSTS, 1);
  enthalpy->set_attrs("model_state",
                      "ice enthalpy (includes sensible heat, latent heat, pressure)",
                      "J kg-1", "");
  variables.add(*enthalpy);
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscLogStage cold, hot;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  ierr = PetscLogStageRegister("Cold", &cold); CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Hot ", &hot); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    Context::Ptr ctx = context_from_options(com, "blatter_test");
    Config::Ptr config = ctx->config();

    ctx->log()->message(2,
                        "BLATTER_TEST: testing the Blatter stress balance solver.\n");

    bool
      usage = options::Bool("-usage", "print usage and stop"),
      help  = options::Bool("-help", "print usage and stop");

    if (usage or help) {
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

    bool compare = options::Bool("-compare", "Compare \"cold\" and \"hot\" runs.");

    ierr = get_grid_from_file(input_file, grid); CHKERRQ(ierr);

    Vars variables;
    allocate_variables(grid, variables);

    EnthalpyConverter EC(config);

    POConstant ocean(grid, config);

    read_input_data(input_file, variables, EC);

    IceModelVec2S melange_back_pressure;
    melange_back_pressure.create(grid, "melange_back_pressure", WITHOUT_GHOSTS);
    melange_back_pressure.set_attrs("boundary_condition",
                                    "melange back pressure fraction",
                                    "1", "");
    melange_back_pressure.set(0.0);

    PetscLogStagePush(cold);
    BlatterStressBalance blatter(grid, EC, config);
    // Initialize the Blatter solver:
    ierr = blatter.init(variables); CHKERRQ(ierr);
    ierr = blatter.update(false, melange_back_pressure); CHKERRQ(ierr);
    PetscLogStagePop();

    if (compare_cold_and_hot) {
      PetscLogStagePush(hot);
      ierr = blatter.update(false, melange_back_pressure); CHKERRQ(ierr);
      PetscLogStagePop();
    }

    // Write results to an output file:
    PIO pio(grid, grid.config.get_string("output_format"));

    ierr = pio.open(output_file, PISM_READWRITE_MOVE); CHKERRQ(ierr);
    ierr = pio.def_time(config.get_string("time_dimension_name"),
                        grid.time->calendar(),
                        grid.time->CF_units_string()); CHKERRQ(ierr);
    ierr = pio.append_time(config.get_string("time_dimension_name"), 0.0);

    std::set<std::string> blatter_vars;
    blatter_vars.insert("u_sigma");
    blatter_vars.insert("v_sigma");
    ierr = blatter.define_variables(blatter_vars, pio, PISM_DOUBLE); CHKERRQ(ierr);
    ierr = blatter.write_variables(blatter_vars, pio); CHKERRQ(ierr);

    ierr = pio.close(); CHKERRQ(ierr);

    ierr = write_data(output_file, variables); CHKERRQ(ierr);

    IceModelVec3 *u, *v;
    IceModelVec2V *velbar;
    ierr = blatter.get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);
    ierr = blatter.get_2D_advective_velocity(velbar); CHKERRQ(ierr);

    ierr = u->write(output_file); CHKERRQ(ierr);
    ierr = v->write(output_file); CHKERRQ(ierr);
    ierr = velbar->write(output_file); CHKERRQ(ierr);


    ierr =  deallocate_variables(variables); CHKERRQ(ierr);

  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
