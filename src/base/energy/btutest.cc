// Copyright (C) 2009, 2010, 2011 Ed Bueler and Constantine Khroulev
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
  "Driver for testing BedrockThermalUnit (without IceModel).\n";

#include "pism_const.hh"
#include "grid.hh"
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "bedrockThermalUnit.hh"


static PetscErrorCode setupIceGridFromFile(string filename, IceGrid &grid) {
  PetscErrorCode ierr;

  PISMIO nc(&grid);
  ierr = nc.get_grid(filename.c_str()); CHKERRQ(ierr);
  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.createDA(); CHKERRQ(ierr);  
  return 0;
}


static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {
  
  PetscErrorCode ierr;
  IceModelVec2S *mask;
  IceModelVec3  *enthalpy;

  mask     = new IceModelVec2S;
  enthalpy = new IceModelVec3;
  
  ierr = mask->create(grid, "mask", true); CHKERRQ(ierr);
  ierr = mask->set_attrs("", "grounded_dragging_floating integer mask",
			      "", ""); CHKERRQ(ierr);
  ierr = variables.add(*mask); CHKERRQ(ierr);

  ierr = enthalpy->create(grid, "enthalpy", false); CHKERRQ(ierr);
  ierr = enthalpy->set_attrs(
     "model_state",
     "ice enthalpy (includes sensible heat, latent heat, pressure)",
     "J kg-1", ""); CHKERRQ(ierr);
  ierr = variables.add(*enthalpy); CHKERRQ(ierr);

  return 0;
}


static PetscErrorCode readIceInfoFromFile(const char *filename, int start,
                                          PISMVars &variables) {
  PetscErrorCode ierr;
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);
    ierr = var->read(filename, start); CHKERRQ(ierr);
    i++;
  }

  return 0;
}


static PetscErrorCode doneWithIceInfo(PISMVars &variables) {
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

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    NCConfigVariable config, overrides;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    // check required options
    vector<string> required;
    required.push_back("-i");
    required.push_back("-o");
    required.push_back("-ys");
    required.push_back("-dt");
    ierr = show_usage_check_req_opts(com, "pclimate", required,
      "  btu_test -i IN.nc -o OUT.nc -ys A -dt B\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "  -ys            start time A (= float) in years\n"
      "  -dt            time step B (= positive float) in years\n"
      ); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    bool override_used;
    ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);

    IceGrid grid(com, rank, size, config);
    
    bool flag;
    PetscReal ys = 0.0, dt_years = 0.0;
    string inname, outname;
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-i", "Input file name", inname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);

      ierr = PISMOptionsReal("-ys", "Start year", ys, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // initialize the computational grid:
    ierr = verbPrintf(2,com, 
	"  initializing grid from NetCDF file %s...\n", inname.c_str()); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);

    grid.year = ys;		
    grid.start_year = ys;
    grid.end_year = ys + dt_years;

    EnthalpyConverter EC(config);

    // allocate IceModelVecs needed by boundary models and put them in a dictionary:
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // read data from a PISM input file
    NCTool nc(grid.com, grid.rank);
    int last_record;
    ierr = nc.open_for_reading(inname.c_str()); CHKERRQ(ierr);
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    last_record -= 1;

    ierr = verbPrintf(2,com, 
        "  reading fields mask,enthalpy from NetCDF file %s to fill fields in PISMVars ...\n",
	inname.c_str()); CHKERRQ(ierr);
    ierr = readIceInfoFromFile(inname.c_str(), last_record, variables); CHKERRQ(ierr);

    // Initialize BTU object:
    BedrockThermalUnit btu(grid, EC, config);

    ierr = btu.update(ys, dt_years); CHKERRQ(ierr);

    ierr = btu.write_model_state(ys, dt_years, outname); CHKERRQ(ierr);

    if (override_used) {
      ierr = verbPrintf(3, com,
        "  recording config overrides in NetCDF file '%s' ...\n",
	outname.c_str()); CHKERRQ(ierr);
      overrides.update_from(config);
      ierr = overrides.write(outname.c_str()); CHKERRQ(ierr);
    }

    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done.\n"); CHKERRQ(ierr);
    
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


