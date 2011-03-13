// Copyright (C) 2011 Ed Bueler and Constantine Khroulev
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
  "Driver for testing PISMBedThermalUnit using Test K.  Sans IceModel.\n"
  "Requires an input file only to get 2D grid information, and otherwise ignors\n"
  "3D and thermodynamic information in the file.\n";

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
  
  // these vertical grid settings are about *input file only* and will not
  //   affect the computation or the output file grid
  grid.compute_vertical_levels();
  grid.printInfo(1);
  grid.printVertLevels(1);

  ierr = grid.createDA(); CHKERRQ(ierr);
  return 0;
}


static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {

  PetscErrorCode ierr;
  IceModelVec2Mask *mask = new IceModelVec2Mask;
  IceModelVec2S *thk = new IceModelVec2S,
                *ghf = new IceModelVec2S;
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

  ierr = ghf->create(grid, "bheatflx", false); CHKERRQ(ierr);
  ierr = ghf->set_attrs("",
                       "upward geothermal flux at bedrock base",
		       "W m-2", ""); CHKERRQ(ierr);
  ierr = ghf->set_glaciological_units("mW m-2");
  ierr = variables.add(*ghf); CHKERRQ(ierr);

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
    required.push_back("-ye");
    required.push_back("-dt");
    required.push_back("-Mz");
    required.push_back("-Mbz");
    required.push_back("-Lbz");
    ierr = show_usage_check_req_opts(com, "btutest", required,
      "  btutest -i IN.nc -o OUT.nc -ys A -ye B -dt C -Mz MM -Mbz NN -Lbz 1000.0\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "  -ys            start year in using Test K\n"
      "  -ye            end year in using Test K\n"
      "  -dt            time step B (= positive float) in years\n"
      "  -Mz            number of ice levels to use; note -Lz option not needed here\n"
      "  -Mbz           number of bedrock thermal layer levels to use\n"
      "  -Lbz 1000.0    depth of bedrock thermal layer (required; Lbz=1000.0 m in Test K)\n"
      ); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid grid(com, rank, size, config);

    bool flag;
    PetscReal dt_years = 0.0, ys, ye;
    PetscInt  Mz;
    string inname, outname;
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-i", "Input file name", inname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ys", "starting time in years", ys, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ye", "starting time in years", ye, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-Mz", "number of ice layers", Mz, flag); CHKERRQ(ierr);
      // -Mbz is checked by IceModelVec3BTU (FIXME)
      // -Lbz is checked by IceModelVec3BTU (FIXME)
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // initialize the computational grid:
    ierr = verbPrintf(2,com,
	"  initializing grid (2D) from NetCDF file %s...\n", inname.c_str()); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);
    grid.start_year = ys;
    grid.end_year   = ye;

    // allocate tools and IceModelVecs
    EnthalpyConverter EC(config);
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

// FIXME  use exactTestK() results and EC to generate needed values in variables

    IceModelVec2Mask *mask;
    mask = dynamic_cast<IceModelVec2Mask*>(variables.get("mask"));
    if (mask == NULL) SETERRQ(2, "mask is not available");
    ierr = mask->set(MASK_SHEET); CHKERRQ(ierr);

    IceModelVec2S *thk;
    thk = dynamic_cast<IceModelVec2S*>(variables.get("thk"));
    if (thk == NULL) SETERRQ(2, "thk is not available");
    ierr = thk->set(3000.0); CHKERRQ(ierr);  // see Test K

    // Initialize BTU object:
    PISMBedThermalUnit btu(grid, EC, config);

    ierr = btu.init(variables); CHKERRQ(ierr);

    // worry about time step
    int  N = (int)ceil((ye - ys) / dt_years);
    PetscReal old_dt_years = dt_years;
    dt_years = (ye - ys) / (double)N;
    ierr = verbPrintf(2,com,
        "  user set timestep of %.4f years ...\n"
        "  reset to %.4f years to get integer number of steps ... \n",
	old_dt_years,dt_years); CHKERRQ(ierr);
    PetscReal max_dt_years;
    ierr = btu.max_timestep(0.0, max_dt_years); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  PISMBedThermalUnit reports max timestep of %.4f years ...\n",
	max_dt_years); CHKERRQ(ierr);

    // actually ask btu to do its time-step
    for (PetscInt n = 0; n < N; n++) {
      // FIXME:  exactTestK() should be queried for the top-of-bedrock temp
      const PetscReal y = ys + dt_years * (double)n;
      ierr = btu.update(y, dt_years); CHKERRQ(ierr);
    }

    IceModelVec2S *ghf;
    ghf = dynamic_cast<IceModelVec2S*>(variables.get("bheatflx"));
    if (ghf == NULL) SETERRQ(1, "bheatflx is not available");

    // get the final G_0 geothermal flux
    ierr = btu.get_upward_geothermal_flux(*ghf); CHKERRQ(ierr);

// FIXME:  compare to the right answer from Test K

    set<string> vars;
    btu.add_vars_to_output("big", vars); // "write everything you can"

    PISMIO pio(&grid);

    ierr = pio.open_for_writing(outname, false, true); CHKERRQ(ierr);
    // append == true and check_dims == true
    ierr = pio.append_time(grid.end_year); CHKERRQ(ierr);
    ierr = btu.define_variables(vars, pio, NC_DOUBLE); CHKERRQ(ierr);
    ierr = pio.close(); CHKERRQ(ierr);

    ierr = btu.write_variables(vars, outname); CHKERRQ(ierr);
    ierr = ghf->write(outname.c_str()); CHKERRQ(ierr);

    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done.\n"); CHKERRQ(ierr);

  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


