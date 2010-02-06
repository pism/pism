// Copyright (C) 2009, 2010 Ed Bueler and Constantine Khroulev
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
  "Driver for testing PISM's boundary models without IceModel.\n";

#include <set>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <petscda.h>
#include "base/pism_const.hh"
#include "base/grid.hh"
#include "base/LocalInterpCtx.hh"
#include "base/nc_util.hh"
#include "base/NCVariable.hh"

#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"
#include "eismint/pgrn_atmosphere.hh"

//! Replaces PetscOptionsEList.
static PetscErrorCode PISMOptionsList(MPI_Comm com, string opt, string description, set<string> choices,
				      string default_value, string &result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  string list, descr;
  PetscTruth opt_set = PETSC_FALSE;

  if (choices.empty()) {
    SETERRQ(1, "PISMOptionsList: empty choices argument");
  }

  list = "[" + *(choices.begin());
  set<string>::iterator j = choices.begin(); j++;
  while (j != choices.end()) {
    list += ", " + (*j++);
  }
  list += "]";

  descr = description + " Choose one of " + list;

  ierr = PetscOptionsString(opt.c_str(), descr.c_str(), "", default_value.c_str(),
			    tmp, TEMPORARY_STRING_LENGTH, &opt_set); CHKERRQ(ierr);

  // return the default value if the option was not set
  if (!opt_set) {
    flag = false;
    result = default_value;
    return 0;
  }

  // return the choice if it is valid and stop if it is not
  if (choices.find(tmp) != choices.end()) {
    flag = true;
    result = tmp;
  } else {
    ierr = PetscPrintf(com, "ERROR: invalid %s argument: \"%s\". Please choose one of %s.\n",
		       opt.c_str(), tmp, list.c_str()); CHKERRQ(ierr);
    ierr = PetscEnd(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMOptionsStrings(string opt, string text, string default_value,
				  vector<string>& result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  PetscTruth opt_set = PETSC_FALSE;

  ierr = PetscOptionsString(opt.c_str(), text.c_str(), "", default_value.c_str(),
			    tmp, TEMPORARY_STRING_LENGTH, &opt_set); CHKERRQ(ierr);

  result.clear();

  if (opt_set) {
    istringstream arg(tmp);
    string word;
    while (getline(arg, word, ','))
      result.push_back(word);

    flag = true;
  } else {
    result.push_back(default_value);
    flag = false;
  }

  return 0;
}
/*
static PetscErrorCode init_atmosphere(IceGrid &grid, const NCConfigVariable &conf, PISMVars &vars,
				      PISMAtmosphereModel* &atmosphere) {
  PetscErrorCode ierr;
  bool opt_set = false;
  string choice;
  PISMAtmosphereModel *tmp = NULL;

  ierr = PetscOptionsBegin(grid.com, "", "PISM options controlling surface, atmosphere and ocean models", ""); CHKERRQ(ierr);

  // Atmosphere models:
  set<string> atmosphere_models;
  atmosphere_models.insert("constant");
  atmosphere_models.insert("greenland");
  atmosphere_models.insert("eismint_greenland");

  ierr = PISMOptionsList(grid.com,
			 "-atmosphere",
			 "Specifies an atmosphere model.",
			 atmosphere_models,
			 "constant",
			 choice,
			 opt_set); CHKERRQ(ierr);

  if        (choice == "constant") {
    atmosphere = new PAConstant(grid, conf, vars);
  } else if (choice == "greenland") {
    atmosphere = new PAFausto(grid, conf, vars);
  } else if (choice == "eismint_greenland") {
    atmosphere = new PA_EISMINT_Greenland(grid, conf, vars);
  } else {
    SETERRQ(1, "Invalid atmosphere model choice (should never happen)");
  }

  return 0;
}
*/

//! Processes command-line options and creates 
static PetscErrorCode init_boundary_models(IceGrid &grid, const NCConfigVariable &conf, PISMVars &vars,
					   PISMSurfaceModel* &surface_model,
					   PISMOceanModel* &ocean_model) {
  PetscErrorCode ierr;
  bool opt_set = false;
  string choice;
  PISMAtmosphereModel *atmosphere = NULL;

  ierr = PetscOptionsBegin(grid.com, "", "PISM options controlling surface, atmosphere and ocean models", ""); CHKERRQ(ierr);

  // Atmosphere models:
  set<string> atmosphere_models;
  atmosphere_models.insert("constant");
  atmosphere_models.insert("greenland");
  atmosphere_models.insert("eismint_greenland");

  ierr = PISMOptionsList(grid.com,
			 "-atmosphere",
			 "Specifies an atmosphere model.",
			 atmosphere_models,
			 "constant",
			 choice,
			 opt_set); CHKERRQ(ierr);

  if        (choice == "constant") {
    atmosphere = new PAConstant(grid, conf, vars);
  } else if (choice == "greenland") {
    atmosphere = new PAFausto(grid, conf, vars);
  } else if (choice == "eismint_greenland") {
    atmosphere = new PA_EISMINT_Greenland(grid, conf, vars);
  } else {
    SETERRQ(1, "Invalid atmosphere model choice (should never happen)");
  }

  // Atmosphere forcing:
  PetscTruth dTforcing_set, temp_ma_anomaly_set, snowprecip_anomaly_set;
  ierr = check_option("-dTforcing", dTforcing_set); CHKERRQ(ierr);
  ierr = check_option("-temp_ma_anomaly", temp_ma_anomaly_set); CHKERRQ(ierr);
  ierr = check_option("-snowprecip_anomaly", snowprecip_anomaly_set); CHKERRQ(ierr);

  if (dTforcing_set || temp_ma_anomaly_set || snowprecip_anomaly_set) {
    PAForcing *forcing = new PAForcing(grid, conf, vars);

    forcing->attach_input(atmosphere);

    atmosphere = forcing;
  }
  
  // Surface models:
  set<string> surface_models;
  surface_models.insert("simple");
  surface_models.insert("constant");
  surface_models.insert("pdd");

  ierr = PISMOptionsList(grid.com,
			 "-surface",
			 "Specifies a surface (snow process, melt) model.",
			 surface_models,
			 "simple",
			 choice,
			 opt_set); CHKERRQ(ierr);

  if        (choice == "simple") {
    surface_model = new PSSimple(grid, conf, vars);
  } else if (choice == "constant") {
    surface_model = new PSConstant(grid, conf, vars);
  } else if (choice == "pdd") {
    surface_model = new PSLocalMassBalance(grid, conf, vars);
  } else {
    SETERRQ(1, "Invalid surface model choice (should never happen)");
  }

  // Surface model modifiers:

  PetscTruth ftt_set;
  ierr = check_option("-force_to_thk", ftt_set); CHKERRQ(ierr);
  if (ftt_set) {
    PSForceThickness *forcing = new PSForceThickness(grid, conf, vars);

    forcing->attach_input(surface_model);

    surface_model = forcing;
  }

  // Ocean models: no choices (yet), so we just use the default (constant) model.

  ocean_model = new POConstant(grid, conf, vars);

  PetscTruth sea_level_forcing = PETSC_FALSE;
  ierr = check_option("-dSLforcing", sea_level_forcing); CHKERRQ(ierr);

  if (sea_level_forcing) {
    POForcing *forcing = new POForcing(grid, conf, vars);

    forcing->attach_input(ocean_model);

    ocean_model = forcing;
  }

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  surface_model->attach_atmosphere_model(atmosphere);

  ierr = surface_model->init(); CHKERRQ(ierr);

  ierr = ocean_model->init(); CHKERRQ(ierr);

  return 0;
}

static PetscErrorCode setupIceGridFromFile(const char *filename, IceGrid &grid) {
  PetscErrorCode ierr;

  NCTool nc(&grid);
  ierr = nc.get_grid(filename); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);  
  return 0;
}

static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {
  
  PetscErrorCode ierr;
  IceModelVec2 *lat, *lon, *mask, *thk, *surfelev, *topg, *acab, *artm, *shelfbasetemp, *shelfbasemassflux;

  lat      = new IceModelVec2;
  lon      = new IceModelVec2;
  mask     = new IceModelVec2;
  thk      = new IceModelVec2;
  surfelev = new IceModelVec2;
  topg     = new IceModelVec2;
  acab     = new IceModelVec2;
  artm     = new IceModelVec2;
  shelfbasetemp     = new IceModelVec2;
  shelfbasemassflux = new IceModelVec2;

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
			 "ice temperature (at the ice surface)",
			 "K",
			 ""); CHKERRQ(ierr);
  ierr = variables.add(*artm); CHKERRQ(ierr);

  ierr = acab->create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab->set_attrs("climate_state", 
			 "ice-equivalent accumulation/ablation rate",
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

static PetscErrorCode readIceInfoFromFile(char *filename, int start,
                                          PISMVars &variables) {
  PetscErrorCode ierr;

  ierr = variables.get("lat")->read(filename, start); CHKERRQ(ierr);
  ierr = variables.get("lon")->read(filename, start); CHKERRQ(ierr);
  ierr = variables.get("mask")->read(filename, start); CHKERRQ(ierr);
  ierr = variables.get("thk")->read(filename, start); CHKERRQ(ierr);
  ierr = variables.get("usurf")->read(filename, start); CHKERRQ(ierr);
  ierr = variables.get("topg")->read(filename, start); CHKERRQ(ierr);

  return 0;
}


static PetscErrorCode doneWithIceInfo(PISMVars &variables) {

  delete variables.get("lat");
  delete variables.get("lon");
  delete variables.get("mask");
  delete variables.get("thk");
  delete variables.get("usurf");
  delete variables.get("topg");

  delete variables.get("artm");
  delete variables.get("acab");

  delete variables.get("shelfbasetemp");
  delete variables.get("shelfbasemassflux");

  return 0;
}


static PetscErrorCode writePCCStateAtTimes(PISMVars &variables,
					   PISMSurfaceModel *surface,
					   PISMOceanModel* ocean,
					   const char *filename, IceGrid* grid,
					   int argc, char *argv[],
					   PetscReal ys, PetscReal ye, PetscReal dt_years,
					   NCConfigVariable &mapping) {

  MPI_Comm com = grid->com;
  PetscErrorCode ierr;
  NCTool nc(grid);
  NCGlobalAttributes global_attrs;
  IceModelVec2 *artm, *acab, *shelfbasetemp, *shelfbasemassflux;

  artm = dynamic_cast<IceModelVec2*>(variables.get("artm"));
  if (artm == NULL) { SETERRQ(1, "artm is not available"); }

  acab = dynamic_cast<IceModelVec2*>(variables.get("acab"));
  if (acab == NULL) { SETERRQ(1, "acab is not available"); }

  shelfbasetemp = dynamic_cast<IceModelVec2*>(variables.get("shelfbasetemp"));
  if (shelfbasetemp == NULL) { SETERRQ(1, "shelfbasetemp is not available"); }

  shelfbasemassflux = dynamic_cast<IceModelVec2*>(variables.get("shelfbasemassflux"));
  if (shelfbasemassflux == NULL) { SETERRQ(1, "shelfbasemassflux is not available"); }

  global_attrs.init("global_attributes", com, grid->rank);
  global_attrs.set_string("Conventions", "CF-1.4");
  global_attrs.set_string("source", string("pcctest3 ") + PISM_Revision);

  // Create a string with space-separated command-line arguments:
  string cmdstr;
  for (int j = 0; j < argc; j++)
    cmdstr += string(" ") + argv[j];
  cmdstr += "\n";

  string history = username_prefix() + cmdstr;

  global_attrs.prepend_history(history);

  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attrs.write(filename); CHKERRQ(ierr);

  PetscInt NN;  // get number of times at which PISM boundary model state is written
  NN = (int) ceil((ye - ys) / dt_years);
  if (NN > 1000)
    SETERRQ(2,"PCCTEST3 ERROR: refuse to write more than 1000 times!");
  if (NN > 50) {
    ierr = PetscPrintf(com,
        "\nPCCTEST3 ATTENTION: writing more than 50 times to '%s'!!\n\n",
        filename); CHKERRQ(ierr);
  }

  PetscScalar use_dt_years = dt_years;
  /*
  PISMGreenlandAtmosCoupler* pdd_pcc = dynamic_cast<PISMGreenlandAtmosCoupler*>(pcc);
  if ((pdd_pcc != NULL) && (dt_years > 1.0)) {
    ierr = verbPrintf(1,com,
        "PCCTEST ATTENTION: Chosen PISM surface model will be asked for results\n"
        "  from one year periods at the start of each desired time subinterval;\n"
        "  full subinterval evaluation is too slow ...\n\n");
      CHKERRQ(ierr);
    use_dt_years = 1.0;
  }
  */

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

    ierr = surface->ice_surface_mass_flux(pccyear, dt_update_years, *acab); CHKERRQ(ierr);
    ierr = acab->write(filename, NC_FLOAT); CHKERRQ(ierr);

    ierr = surface->ice_surface_temperature(pccyear, dt_update_years, *artm); CHKERRQ(ierr);
    ierr = artm->write(filename, NC_FLOAT); CHKERRQ(ierr);

    ierr = surface->write_diagnostic_fields(pccyear, dt_update_years, filename); CHKERRQ(ierr);

    ierr = ocean->shelf_base_temperature(pccyear, dt_update_years, *shelfbasetemp); CHKERRQ(ierr);
    ierr = shelfbasetemp->write(filename, NC_FLOAT); CHKERRQ(ierr);

    ierr = ocean->shelf_base_mass_flux(pccyear, dt_update_years, *shelfbasemassflux); CHKERRQ(ierr);
    ierr = shelfbasemassflux->write(filename, NC_FLOAT); CHKERRQ(ierr);

    ierr = ocean->write_diagnostic_fields(pccyear, dt_update_years, filename); CHKERRQ(ierr);
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
    IceGrid grid(com, rank, size);
    char inname[PETSC_MAX_PATH_LEN], outname[PETSC_MAX_PATH_LEN];

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    // check required options
    vector<string> required;
    required.push_back("-i");
    required.push_back("-o");
    required.push_back("-ys");
    required.push_back("-ye");
    required.push_back("-dt");
    ierr = show_usage_check_req_opts(com, "pcctest", required,
      "  pcctest -i IN.nc -o OUT.nc -ys A -ye B -dt C [-atmosphere <name> -surface <name>] [OTHER PISM & PETSc OPTIONS]\n\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "  -ys            start time A (= float) in years\n"
      "  -ye            end time B (= float), B > A, in years\n"
      "  -dt            time step C (= positive float) in years\n"
      "and set up the models:\n"
      "  -atmosphere    Chooses an atmosphere model (one of [constant, greenland])\n\n"
      "  -surface       Chooses a surface model (one of [simple, pdd]) \n"
      "  -pdd_greenland Sets PDD parameters using formulas (6) and (7) in [Faustoetal2009]\n"
      ); CHKERRQ(ierr);

    ierr = verbPrintf(2,
      com,"PCCTEST %s (test of PISM boundary models offline from IceModel)\n",
      PISM_Revision); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    // initialize the computational grid:
    ierr = PetscOptionsGetString(PETSC_NULL, "-i", inname, 
                                 PETSC_MAX_PATH_LEN, NULL); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, 
             "  initializing grid from NetCDF file %s...\n", inname); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);

    mapping.init("mapping", com, rank);

    // Process -ys, -ye, -dt *before* initializing boundary models
    PetscReal ys = 0.0, ye = 0.0, dt_years = 0.0;
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ys", &ys, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ye", &ye, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-dt", &dt_years, NULL); CHKERRQ(ierr);
    // These values may be used by surface, atmosphere of ocean models:
    grid.year = ys;		
    grid.start_year = ys;
    grid.end_year = ye;

    // allocate IceModelVecs needed by boundary models and put them in a dictionary:
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // read data from a PISM input file (including the projection parameters)
    NCTool nc(&grid);
    int last_record;
    bool mapping_exists;
    ierr = nc.open_for_reading(inname); CHKERRQ(ierr);
    ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    if (mapping_exists) {
      ierr = mapping.read(inname); CHKERRQ(ierr);
      ierr = mapping.print(); CHKERRQ(ierr);
    }
    last_record -= 1;

    ierr = verbPrintf(2,com, 
             "  reading fields lat,lon,mask,thk,topg,usurf from NetCDF file %s\n"
             "    to fill fields in PISMVars ...\n",
             inname); CHKERRQ(ierr);

    ierr = readIceInfoFromFile(inname, last_record, variables); CHKERRQ(ierr);

    // Initialize boundary models:
    PISMSurfaceModel *surface = NULL;
    PISMOceanModel   *ocean   = NULL;
    ierr = init_boundary_models(grid, config, variables, surface, ocean); CHKERRQ(ierr);
    
    ierr = PetscOptionsGetString(PETSC_NULL, "-o", outname, 
                                 PETSC_MAX_PATH_LEN, NULL); CHKERRQ(ierr);

    ierr = verbPrintf(2,
      com, "  writing boundary model states to NetCDF file '%s'...\n",
      outname); CHKERRQ(ierr);

    ierr = writePCCStateAtTimes(variables, surface, ocean, outname, &grid, argc, argv,
				ys, ye, dt_years,
                                mapping); CHKERRQ(ierr);

    delete surface;
    delete ocean;
    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);
    
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


