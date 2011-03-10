// Copyright (C) 2011 Ed Bueler
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
  "Driver for testing PISMBedThermalUnit (without IceModel).\n";

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
  IceModelVec2S *mask = new IceModelVec2S,
                *thk = new IceModelVec2S,
                *ghf_result = new IceModelVec2S;
  IceModelVec3  *enthalpy = new IceModelVec3;
  
  ierr = mask->create(grid, "mask", true); CHKERRQ(ierr);
  ierr = mask->set_attrs("", "grounded_dragging_floating integer mask",
			      "", ""); CHKERRQ(ierr);
  ierr = variables.add(*mask); CHKERRQ(ierr);

  ierr = thk->create(grid, "thk", true); CHKERRQ(ierr);
  ierr = thk->set_attrs("", "ice thickness",
			      "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = variables.add(*thk); CHKERRQ(ierr);

  ierr = enthalpy->create(grid, "enthalpy", false); CHKERRQ(ierr);
  ierr = enthalpy->set_attrs(
     "",
     "ice enthalpy (includes sensible heat, latent heat, pressure)",
     "J kg-1", ""); CHKERRQ(ierr);
  ierr = variables.add(*enthalpy); CHKERRQ(ierr);

  ierr = ghf_result->create(grid, "bheatflx_at_ice_base", false); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_upward_heat_flux ?
  ierr = ghf_result->set_attrs("",
                       "upward geothermal flux at bedrock surface, at ice base",
		       "W m-2", ""); CHKERRQ(ierr);
  ierr = ghf_result->set_glaciological_units("mW m-2");
  ierr = variables.add(*ghf_result); CHKERRQ(ierr);

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


static PetscErrorCode writeState(PISMVars &variables,
                                 const char *filename, IceGrid* grid) {

  PetscErrorCode ierr;

  MPI_Comm com = grid->com;
  NCGlobalAttributes global_attrs;
  global_attrs.init("global_attributes", com, grid->rank);
  global_attrs.set_string("Conventions", "CF-1.4");
  global_attrs.set_string("source", string("btutest ") + PISM_Revision);
  // Create a string with space-separated command-line arguments:
  string history = pism_username_prefix() + pism_args_string();
  global_attrs.prepend_history(history);

  PISMIO nc(grid);
  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = global_attrs.write(filename); CHKERRQ(ierr);

  ierr = verbPrintf(2,com,"\n  writing result to %s ..", filename); CHKERRQ(ierr);

  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);
    ierr = var->write(filename, NC_FLOAT); CHKERRQ(ierr);
    i++;
  }

  ierr = verbPrintf(2,com,"\n"); CHKERRQ(ierr);

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
    required.push_back("-dt");
    ierr = show_usage_check_req_opts(com, "btutest", required,
      "  btutest -i IN.nc -o OUT.nc -ys A -dt B\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "  -dt            time step B (= positive float) in years\n"
      ); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    bool override_used;
    ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);

    IceGrid grid(com, rank, size, config);
    
    bool flag;
    PetscReal dt_years = 0.0;
    string inname, outname;
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-i", "Input file name", inname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // initialize the computational grid:
    ierr = verbPrintf(2,com, 
	"  initializing grid from NetCDF file %s...\n", inname.c_str()); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);
    grid.end_year = grid.start_year + dt_years;

    // allocate tools and IceModelVecs
    EnthalpyConverter EC(config);
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // read data from a PISM input file
    NCTool nc(grid.com, grid.rank);
    int last_record;
    ierr = nc.open_for_reading(inname.c_str()); CHKERRQ(ierr);
    ierr = nc.get_nrecords(last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    last_record -= 1;
    ierr = verbPrintf(2,com, 
        "  reading fields mask,thk,enthalpy from NetCDF file %s ...\n",
	inname.c_str()); CHKERRQ(ierr);
    ierr = readIceInfoFromFile(inname.c_str(), last_record, variables); CHKERRQ(ierr);

    // Initialize BTU object:
    PISMBedThermalUnit btu(grid, EC, config);

    ierr = verbPrintf(2,com, 
        "  user set timestep of %.4f years ...\n",
	dt_years); CHKERRQ(ierr);
    PetscReal max_dt_years;
    ierr = btu.max_timestep(0.0, max_dt_years); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, 
        "  PISMBedThermalUnit reports max timestep of %.4f years ...\n",
	max_dt_years); CHKERRQ(ierr);

    // FIXME:  use btu.get_upward_geothermal_flux(); requires actually having the
    //         IceModelVec2S for the result; use variables.get("???") ?
    
    ierr = writeState(variables, outname.c_str(), &grid); CHKERRQ(ierr);

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


