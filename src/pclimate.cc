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
  "Driver for testing PISM's boundary (surface and shelf-base) models without IceModel.\n";

#include <set>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <petscda.h>
#include "pism_const.hh"
#include "grid.hh"
#include "LocalInterpCtx.hh"
#include "PISMIO.hh"
#include "NCVariable.hh"

#include "coupler/PCFactory.hh"
#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"
#include "eismint/pgrn_atmosphere.hh"

static void create_pa_eismint_greenland(IceGrid& g, const NCConfigVariable& conf,
					PISMAtmosphereModel* &result) {
  result = new PA_EISMINT_Greenland(g, conf);
}


static PetscErrorCode setupIceGridFromFile(string filename, IceGrid &grid) {
  PetscErrorCode ierr;

  PISMIO nc(&grid);
  ierr = nc.get_grid("land_ice_thickness", filename.c_str()); CHKERRQ(ierr);
  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.createDA(); CHKERRQ(ierr);  
  return 0;
}

static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {
  
  PetscErrorCode ierr;
  IceModelVec2S *lat, *lon, *mask, *thk, *surfelev, *topg, *acab, *artm, *shelfbasetemp, *shelfbasemassflux;

  lat      = new IceModelVec2S;
  lon      = new IceModelVec2S;
  mask     = new IceModelVec2S;
  thk      = new IceModelVec2S;
  surfelev = new IceModelVec2S;
  topg     = new IceModelVec2S;
  
  // the following are allocated by the pclimate code, but may or may not
  //   actually be read by PISMAtmosphereModel *atmosphere and PISMOceanModel *ocean
  acab     = new IceModelVec2S;
  artm     = new IceModelVec2S;
  shelfbasetemp     = new IceModelVec2S;
  shelfbasemassflux = new IceModelVec2S;

  ierr = lat->create(grid, "lat", true); CHKERRQ(ierr);
  ierr = lat->set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);
  ierr = variables.add(*lat); CHKERRQ(ierr);

  ierr = lon->create(grid, "lon", true); CHKERRQ(ierr);
  ierr = lon->set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);
  ierr = variables.add(*lon); CHKERRQ(ierr);

  ierr = mask->create(grid, "mask", true); CHKERRQ(ierr);
  ierr = mask->set_attrs("", "grounded_dragging_floating integer mask",
			      "", ""); CHKERRQ(ierr);
  ierr = variables.add(*mask); CHKERRQ(ierr);

  ierr = thk->create(grid, "thk", true); CHKERRQ(ierr);
  ierr = thk->set_attrs("", "land ice thickness",
		             "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = variables.add(*thk); CHKERRQ(ierr);

  ierr = surfelev->create(grid, "usurf", true); CHKERRQ(ierr);
  ierr = surfelev->set_attrs("", "ice upper surface elevation",
		                  "m", "surface_altitude"); CHKERRQ(ierr);
  ierr = variables.add(*surfelev); CHKERRQ(ierr);

  ierr = topg->create(grid, "topg", true); CHKERRQ(ierr);
  ierr = topg->set_attrs("", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);
  ierr = variables.add(*topg); CHKERRQ(ierr);

  ierr = artm->create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm->set_attrs("climate_state",
			 "annual average ice surface temperature, below firn processes",
			 "K",
			 ""); CHKERRQ(ierr);
  ierr = variables.add(*artm); CHKERRQ(ierr);

  ierr = acab->create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab->set_attrs("climate_state", 
			 "ice-equivalent surface mass balance (accumulation/ablation) rate",
			 "m s-1", 
			 ""); CHKERRQ(ierr);
  ierr = acab->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab->write_in_glaciological_units = true;
  ierr = variables.add(*acab); CHKERRQ(ierr);

  ierr = shelfbasetemp->create(grid, "shelfbasetemp", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbasetemp->set_attrs(
				 "climate_state", "absolute temperature at ice shelf base",
				 "K", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_temperature
  ierr = variables.add(*shelfbasetemp); CHKERRQ(ierr);

  // ice mass balance rate at the base of the ice shelf; sign convention for vshelfbasemass
  //   matches standard sign convention for basal melt rate of grounded ice
  ierr = shelfbasemassflux->create(grid, "shelfbasemassflux", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbasemassflux->set_attrs("climate_state",
				     "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
				     "m s-1", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  shelfbasemassflux->write_in_glaciological_units = true;
  ierr = shelfbasemassflux->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = variables.add(*shelfbasemassflux); CHKERRQ(ierr);


  return 0;
}

static PetscErrorCode readIceInfoFromFile(const char *filename, int start,
                                          PISMVars &variables) {
  PetscErrorCode ierr;
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  // Remove artm, acab, shelfbasemassflux and shelfbasetemp: they are filled by
  // surface and ocean models and aren't necessarily read from files.
  vars.erase("artm");
  vars.erase("acab");
  vars.erase("shelfbasemassflux");
  vars.erase("shelfbasetemp");

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


static PetscErrorCode writePCCStateAtTimes(PISMVars &variables,
					   PISMSurfaceModel *surface,
					   PISMOceanModel* ocean,
					   const char *filename, IceGrid* grid,
					   PetscReal ys, PetscReal ye, PetscReal dt_years,
					   NCConfigVariable &mapping) {

  MPI_Comm com = grid->com;
  PetscErrorCode ierr;
  PISMIO nc(grid);
  NCGlobalAttributes global_attrs;
  IceModelVec2S *usurf, *artm, *acab, *shelfbasetemp, *shelfbasemassflux;

  usurf = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (usurf == NULL) { SETERRQ(1, "usurf is not available"); }

  artm = dynamic_cast<IceModelVec2S*>(variables.get("artm"));
  if (artm == NULL) { SETERRQ(1, "artm is not available"); }

  acab = dynamic_cast<IceModelVec2S*>(variables.get("acab"));
  if (acab == NULL) { SETERRQ(1, "acab is not available"); }

  shelfbasetemp = dynamic_cast<IceModelVec2S*>(variables.get("shelfbasetemp"));
  if (shelfbasetemp == NULL) { SETERRQ(1, "shelfbasetemp is not available"); }

  shelfbasemassflux = dynamic_cast<IceModelVec2S*>(variables.get("shelfbasemassflux"));
  if (shelfbasemassflux == NULL) { SETERRQ(1, "shelfbasemassflux is not available"); }

  global_attrs.init("global_attributes", com, grid->rank);
  global_attrs.set_string("Conventions", "CF-1.4");
  global_attrs.set_string("source", string("pclimate ") + PISM_Revision);

  // Create a string with space-separated command-line arguments:
  string history = pism_username_prefix() + pism_args_string();

  global_attrs.prepend_history(history);

  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attrs.write(filename); CHKERRQ(ierr);

  PetscInt NN;  // get number of times at which PISM boundary model state is written
  NN = (int) ceil((ye - ys) / dt_years);
  if (NN > 1000)
    SETERRQ(2,"PCLIMATE ERROR: refuse to write more than 1000 times!");
  if (NN > 50) {
    ierr = PetscPrintf(com,
        "\nPCLIMATE ATTENTION: writing more than 50 times to '%s'!!\n\n",
        filename); CHKERRQ(ierr);
  }

  PetscScalar use_dt_years = dt_years;

  set<string> vars_to_write;
  surface->add_vars_to_output("big", vars_to_write);
  ocean->add_vars_to_output("big", vars_to_write);

  // write the states
  for (PetscInt k=0; k < NN; k++) {
    // use original dt_years to get correct subinterval starts:
    const PetscReal pccyear = ys + k * dt_years; 
    ierr = nc.open_for_writing(filename, true, false); CHKERRQ(ierr); // append=true,check_dims=false
    ierr = nc.append_time(pccyear); CHKERRQ(ierr);
    
    PetscScalar dt_update_years = PetscMin(use_dt_years, ye - pccyear);

    char timestr[TEMPORARY_STRING_LENGTH];
    snprintf(timestr, sizeof(timestr), 
        "  boundary models updated for [%11.3f a,%11.3f a] ...", 
        pccyear, pccyear + dt_update_years);
    ierr = verbPrintf(2,com,"."); CHKERRQ(ierr);
    ierr = verbPrintf(3,com,"\n%s writing result to %s ..",timestr,filename); CHKERRQ(ierr);
    strncat(timestr,"\n",1);
    ierr = nc.write_history(timestr); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = usurf->write(filename, NC_FLOAT); CHKERRQ(ierr);

    // update surface and ocean models' outputs:
    ierr = surface->ice_surface_mass_flux(pccyear, dt_update_years, *acab); CHKERRQ(ierr);

    ierr = surface->ice_surface_temperature(pccyear, dt_update_years, *artm); CHKERRQ(ierr);

    ierr = ocean->shelf_base_temperature(pccyear, dt_update_years, *shelfbasetemp); CHKERRQ(ierr);

    ierr = ocean->shelf_base_mass_flux(pccyear, dt_update_years, *shelfbasemassflux); CHKERRQ(ierr);

    // ask ocean and surface models to write variables:
    ierr = surface->write_variables(vars_to_write, filename); CHKERRQ(ierr);
    ierr = ocean->write_variables(vars_to_write, filename); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,com,"\n"); CHKERRQ(ierr);

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
    NCConfigVariable config, overrides, mapping;
    string inname, outname;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com,
      "PCLIMATE %s (surface and shelf-base boundary-models-only mode)\n",
      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    // check required options
    vector<string> required;
    required.push_back("-i");
    required.push_back("-o");
    required.push_back("-ys");
    required.push_back("-ye");
    required.push_back("-dt");
    ierr = show_usage_check_req_opts(com, "pclimate", required,
      "  pclimate -i IN.nc -o OUT.nc -ys A -ye B -dt C [-atmosphere <name> -surface <name>] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "  -ys            start time A (= float) in years\n"
      "  -ye            end time B (= float), B > A, in years\n"
      "  -dt            time step C (= positive float) in years\n"
      "and set up the models:\n"
      "  -atmosphere    Chooses an atmosphere model; see User's Manual\n"
      "  -surface       Chooses a surface model; see User's Manual\n"
      ); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    bool override_used;
    ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);

    // set an un-documented (!) flag to limit time-steps to 1 year.
    config.set_flag("pdd_limit_timestep", true);

    IceGrid grid(com, rank, size, config);
    
    bool flag;
    PetscReal ys = 0.0, ye = 0.0, dt_years = 0.0;
    ierr = PetscOptionsBegin(grid.com, "", "PCLIMATE options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-i", "Input file name",  inname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);

      ierr = PISMOptionsReal("-ys", "Start year", ys, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ye", "End year",   ye, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // initialize the computational grid:
    ierr = verbPrintf(2,com, 
		      "  initializing grid from NetCDF file %s...\n", inname.c_str()); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);

    mapping.init("mapping", com, rank);

    // These values may be used by surface, atmosphere of ocean models:
    grid.year = ys;		
    grid.start_year = ys;
    grid.end_year = ye;

    // allocate IceModelVecs needed by boundary models and put them in a dictionary:
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // read data from a PISM input file (including the projection parameters)
    NCTool nc(grid.com, grid.rank);
    int last_record;
    bool mapping_exists;
    ierr = nc.open_for_reading(inname.c_str()); CHKERRQ(ierr);
    ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
    ierr = nc.get_nrecords(last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    if (mapping_exists) {
      ierr = mapping.read(inname.c_str()); CHKERRQ(ierr);
      ierr = mapping.print(); CHKERRQ(ierr);
    }
    last_record -= 1;

    ierr = verbPrintf(2,com, 
             "  reading fields lat,lon,mask,thk,topg,usurf from NetCDF file %s\n"
             "    to fill fields in PISMVars ...\n",
		      inname.c_str()); CHKERRQ(ierr);

    ierr = readIceInfoFromFile(inname.c_str(), last_record, variables); CHKERRQ(ierr);

    // Initialize boundary models:
    PAFactory pa(grid, config);
    PISMAtmosphereModel *atmosphere;
    pa.add_model("eismint_greenland", &create_pa_eismint_greenland);

    PSFactory ps(grid, config);
    PISMSurfaceModel *surface;

    POFactory po(grid, config);
    PISMOceanModel *ocean;

    ierr = PetscOptionsBegin(grid.com, "", "PISM Boundary Models", ""); CHKERRQ(ierr);

    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);

    surface->attach_atmosphere_model(atmosphere);
    ierr = surface->init(variables); CHKERRQ(ierr);
    ierr = ocean->init(variables); CHKERRQ(ierr);

    ierr = PetscOptionsEnd(); CHKERRQ(ierr);  // done initializing boundary models.

    ierr = verbPrintf(2, com,
        "writing boundary model states to NetCDF file '%s' ...\n",
        outname.c_str()); CHKERRQ(ierr);

    ierr = writePCCStateAtTimes(variables, surface, ocean,
                                outname.c_str(), &grid,
				ys, ye, dt_years,
                                mapping); CHKERRQ(ierr);

    if (override_used) {
      ierr = verbPrintf(3, com,
        "  recording config overrides in NetCDF file '%s' ...\n",
	outname.c_str()); CHKERRQ(ierr);
      overrides.update_from(config);
      ierr = overrides.write(outname.c_str()); CHKERRQ(ierr);
    }

    delete surface;
    delete ocean;
    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "done.\n"); CHKERRQ(ierr);
    
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


