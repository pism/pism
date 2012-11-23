// Copyright (C) 2011, 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "BlatterStressBalance.hh"
#include "IceGrid.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "flowlaws.hh"
#include "enthalpyConverter.hh"
#include "basal_resistance.hh"
#include "pism_options.hh"
#include "PISMTime.hh"

static PetscErrorCode get_grid_from_file(string filename, IceGrid &grid) {
  PetscErrorCode ierr;

  PIO nc(grid, "guess_format");

  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
  ierr = nc.inq_grid("bedrock_altitude", &grid, NOT_PERIODIC); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  grid.compute_nprocs();
  grid.compute_ownership_ranges();

  ierr = grid.createDA(); CHKERRQ(ierr);

  ierr = grid.printInfo(1); CHKERRQ(ierr);

  return 0;
}

//! \brief Read data from an input file.
static PetscErrorCode read_input_data(string filename, PISMVars &variables) {
  PetscErrorCode ierr;
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  for (set<string>::iterator i = vars.begin(); i != vars.end(); ++i) {
    ierr = variables.get(*i)->regrid(filename, true); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Write data to an output file.
static PetscErrorCode write_data(string filename, PISMVars &variables) {
  PetscErrorCode ierr;
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  for (set<string>::iterator i = vars.begin(); i != vars.end(); ++i) {
    ierr = variables.get(*i)->write(filename); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Allocate IceModelVec2S variables.
static PetscErrorCode allocate_variables(IceGrid &grid, PISMVars &variables) {
  PetscErrorCode ierr;
  IceModelVec2S *surfelev, *topg, *tauc;

  surfelev = new IceModelVec2S;
  topg     = new IceModelVec2S;
  tauc     = new IceModelVec2S;

  ierr = surfelev->create(grid, "usurf", false); CHKERRQ(ierr);
  ierr = surfelev->set_attrs("", "ice upper surface elevation",
                             "m", "surface_altitude"); CHKERRQ(ierr);
  ierr = variables.add(*surfelev); CHKERRQ(ierr);

  ierr = topg->create(grid, "topg", false); CHKERRQ(ierr);
  ierr = topg->set_attrs("", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);
  ierr = variables.add(*topg); CHKERRQ(ierr);

  ierr = tauc->create(grid, "tauc", false); CHKERRQ(ierr);
  ierr = tauc->set_attrs("diagnostic",
                         "yield stress for basal till (plastic or pseudo-plastic model)",
                         "Pa", ""); CHKERRQ(ierr);
  ierr = variables.add(*tauc); CHKERRQ(ierr);

  return 0;
}

//! \brief De-allocate IceModelVec2S variables.
static PetscErrorCode deallocate_variables(PISMVars &variables) {
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);
    delete var;
    i++;
  }

  return 0;
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);
    ierr = set_config_from_options(com, config); CHKERRQ(ierr);
    ierr = setVerbosityLevel(2);

    ierr = verbPrintf(2, com,
                      "BLATTER_TEST: testing the Blatter stress balance solver.\n"); CHKERRQ(ierr);

    PetscBool usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of BLATTER_TEST:\n"
                  "  run siafd_test -i input.nc -o output.nc\n"
                  "\n");
    }

    IceGrid grid(com, rank, size, config);

    string input_file, output_file = "blatter_test.nc";
    ierr = PetscOptionsBegin(grid.com, "", "BLATTER_TEST options", ""); CHKERRQ(ierr);
    {
      bool flag;
      ierr = PISMOptionsString("-i", "Set the input file name",
                               input_file, flag); CHKERRQ(ierr);
      if (! flag) {
        PetscPrintf(grid.com, "BLATTER_TEST ERROR: -i is required.\n");
        PISMEnd();
      }
      ierr = PISMOptionsString("-o", "Set the output file name",
                               output_file, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    ierr = get_grid_from_file(input_file, grid); CHKERRQ(ierr);

    PISMVars variables;
    ierr = allocate_variables(grid, variables); CHKERRQ(ierr);

    EnthalpyConverter EC(config);
    ThermoGlenArrIce ice(grid.com, "", config, &EC);

    // This is never used (but it is a required argument of the
    // PISMStressBalance constructor).
    IceBasalResistancePlasticLaw basal(
           config.get("plastic_regularization", "1/year", "1/second"),
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold", "m/year", "m/second"));

    // POConstant ocean(grid, config);

    BlatterStressBalance blatter(grid, NULL, config);

    ierr =  read_input_data(input_file, variables); CHKERRQ(ierr);

    // Initialize the Blatter solver:
    ierr = blatter.init(variables); CHKERRQ(ierr);

    ierr = blatter.update(false); CHKERRQ(ierr);

    // Write results to an output file:
    PIO pio(grid, grid.config.get_string("output_format"));

    ierr = pio.open(output_file, PISM_WRITE); CHKERRQ(ierr);
    ierr = pio.def_time(config.get_string("time_dimension_name"),
                        config.get_string("calendar"),
                        grid.time->CF_units()); CHKERRQ(ierr);
    ierr = pio.append_time(config.get_string("time_dimension_name"), 0.0);
    ierr = pio.close(); CHKERRQ(ierr);

    ierr =  write_data(output_file, variables); CHKERRQ(ierr);

    IceModelVec3 *u, *v, *w;
    ierr =  blatter.get_3d_velocity(u, v, w); CHKERRQ(ierr);

    ierr =  u->write(output_file); CHKERRQ(ierr);
    ierr =  v->write(output_file); CHKERRQ(ierr);
    ierr =  w->write(output_file); CHKERRQ(ierr);

    ierr =  deallocate_variables(variables); CHKERRQ(ierr);

  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
