// Copyright (C) 2009, 2010, 2011, 2012 Ed Bueler and Constantine Khroulev
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
#include <petscdmda.h>
#include "pism_const.hh"
#include "pism_options.hh"
#include "IceGrid.hh"
#include "LocalInterpCtx.hh"
#include "PIO.hh"
#include "NCVariable.hh"
#include "Timeseries.hh"

#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"


static PetscErrorCode setupIceGridFromFile(string filename, IceGrid &grid) {
  PetscErrorCode ierr;

  PIO nc(grid, "guess_format");

  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
  // filename should point to a PISM output file, which is guaranteed to have
  // land_ice_thickness in it.
  ierr = nc.inq_grid("land_ice_thickness", &grid, NOT_PERIODIC); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.createDA(); CHKERRQ(ierr);
  return 0;
}

static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {
  
  PetscErrorCode ierr;
  IceModelVec2S *lat, *lon, *mask, *thk, *surfelev, *topg, *climatic_mass_balance, *ice_surface_temp, *shelfbasetemp, *shelfbasemassflux;

  lat      = new IceModelVec2S;
  lon      = new IceModelVec2S;
  mask     = new IceModelVec2S;
  thk      = new IceModelVec2S;
  surfelev = new IceModelVec2S;
  topg     = new IceModelVec2S;
  
  // the following are allocated by the pclimate code, but may or may not
  //   actually be read by PISMAtmosphereModel *atmosphere and PISMOceanModel *ocean
  climatic_mass_balance     = new IceModelVec2S;
  ice_surface_temp     = new IceModelVec2S;
  shelfbasetemp     = new IceModelVec2S;
  shelfbasemassflux = new IceModelVec2S;

  ierr = lat->create(grid, "lat", true); CHKERRQ(ierr);
  ierr = lat->set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);
  lat->time_independent = true;
  ierr = variables.add(*lat); CHKERRQ(ierr);

  ierr = lon->create(grid, "lon", true); CHKERRQ(ierr);
  ierr = lon->set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);
  lon->time_independent = true;
  ierr = variables.add(*lon); CHKERRQ(ierr);

  ierr = mask->create(grid, "mask", true); CHKERRQ(ierr);
  ierr = mask->set_attrs("", "grounded_dragging_floating integer mask",
			      "", ""); CHKERRQ(ierr);
  mask->time_independent = true;
  ierr = variables.add(*mask); CHKERRQ(ierr);

  ierr = thk->create(grid, "thk", true); CHKERRQ(ierr);
  ierr = thk->set_attrs("", "land ice thickness",
		             "m", "land_ice_thickness"); CHKERRQ(ierr);
  thk->time_independent = true;
  ierr = variables.add(*thk); CHKERRQ(ierr);

  ierr = surfelev->create(grid, "usurf", true); CHKERRQ(ierr);
  ierr = surfelev->set_attrs("", "ice upper surface elevation",
		                  "m", "surface_altitude"); CHKERRQ(ierr);
  surfelev->time_independent = true;
  ierr = variables.add(*surfelev); CHKERRQ(ierr);

  ierr = topg->create(grid, "topg", true); CHKERRQ(ierr);
  ierr = topg->set_attrs("", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);
  topg->time_independent = true;
  ierr = variables.add(*topg); CHKERRQ(ierr);

  ierr = ice_surface_temp->create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
  ierr = ice_surface_temp->set_attrs("climate_state",
                                     "annual average ice surface temperature, below firn processes",
                                     "K",
                                     ""); CHKERRQ(ierr);
  ierr = variables.add(*ice_surface_temp); CHKERRQ(ierr);

  ierr = climatic_mass_balance->create(grid, "climatic_mass_balance", false); CHKERRQ(ierr);
  ierr = climatic_mass_balance->set_attrs("climate_state",
                                          "ice-equivalent surface mass balance (accumulation/ablation) rate",
                                          "m s-1",
                                          ""); CHKERRQ(ierr);
  ierr = climatic_mass_balance->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  climatic_mass_balance->write_in_glaciological_units = true;
  ierr = variables.add(*climatic_mass_balance); CHKERRQ(ierr);

  ierr = shelfbasetemp->create(grid, "shelfbtemp", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbasetemp->set_attrs(
				 "climate_state", "absolute temperature at ice shelf base",
				 "K", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_temperature
  ierr = variables.add(*shelfbasetemp); CHKERRQ(ierr);

  // ice mass balance rate at the base of the ice shelf; sign convention for vshelfbasemass
  //   matches standard sign convention for basal melt rate of grounded ice
  ierr = shelfbasemassflux->create(grid, "shelfbmassflux", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbasemassflux->set_attrs("climate_state",
				     "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
				     "m s-1", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  shelfbasemassflux->write_in_glaciological_units = true;
  ierr = shelfbasemassflux->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = variables.add(*shelfbasemassflux); CHKERRQ(ierr);

  return 0;
}

static PetscErrorCode readIceInfoFromFile(string filename, int start,
                                          PISMVars &variables) {
  PetscErrorCode ierr;
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();

  // Remove ice_surface_temp, climatic_mass_balance, shelfbasemassflux and shelfbasetemp: they are filled by
  // surface and ocean models and aren't necessarily read from files.
  vars.erase("ice_surface_temp");
  vars.erase("climatic_mass_balance");
  vars.erase("shelfbmassflux");
  vars.erase("shelfbtemp");

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
    delete variables.get(*i);
    i++;
  }

  return 0;
}


static PetscErrorCode writePCCStateAtTimes(PISMVars &variables,
					   PISMSurfaceModel *surface,
					   PISMOceanModel* ocean,
					   string filename, IceGrid& grid,
                                           vector<double> times,
					   NCConfigVariable &mapping) {

  MPI_Comm com = grid.com;
  PetscErrorCode ierr;
  PIO nc(grid, grid.config.get_string("output_format"));
  NCGlobalAttributes global_attrs;
  IceModelVec2S *usurf, *ice_surface_temp, *climatic_mass_balance, *shelfbasetemp, *shelfbasemassflux;

  usurf = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (usurf == NULL) { SETERRQ(com, 1, "usurf is not available"); }

  ice_surface_temp = dynamic_cast<IceModelVec2S*>(variables.get("ice_surface_temp"));
  if (ice_surface_temp == NULL) { SETERRQ(com, 1, "ice_surface_temp is not available"); }

  climatic_mass_balance = dynamic_cast<IceModelVec2S*>(variables.get("climatic_mass_balance"));
  if (climatic_mass_balance == NULL) { SETERRQ(com, 1, "climatic_mass_balance is not available"); }

  shelfbasetemp = dynamic_cast<IceModelVec2S*>(variables.get("shelfbtemp"));
  if (shelfbasetemp == NULL) { SETERRQ(com, 1, "shelfbasetemp is not available"); }

  shelfbasemassflux = dynamic_cast<IceModelVec2S*>(variables.get("shelfbmassflux"));
  if (shelfbasemassflux == NULL) { SETERRQ(com, 1, "shelfbasemassflux is not available"); }

  global_attrs.init("global_attributes", com, grid.rank);
  global_attrs.set_string("Conventions", "CF-1.5");
  global_attrs.set_string("source", string("pclimate ") + PISM_Revision);

  // Create a string with space-separated command-line arguments:
  string history = pism_username_prefix(com) + pism_args_string();

  global_attrs.prepend_history(history);

  ierr = nc.open(filename, PISM_WRITE); CHKERRQ(ierr);
  // append == false, check_dims == true

  ierr = mapping.write(nc); CHKERRQ(ierr);
  ierr = global_attrs.write(nc); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  if (times.size() > 1000) {
    PetscPrintf(grid.com, "PCLIMATE ERROR: refuse to write more than 1000 times!");
    PISMEnd();
  }

  if (times.size() > 50) {
    ierr = PetscPrintf(com,
        "\nPCLIMATE ATTENTION: writing more than 50 times to '%s'!!\n\n",
        filename.c_str()); CHKERRQ(ierr);
  }

  DiagnosticTimeseries sea_level(&grid, "sea_level", grid.config.get_string("time_dimension_name"));
  sea_level.set_units("m", "m");
  sea_level.set_dimension_units(grid.time->units(), "");
  sea_level.output_filename = filename;
  sea_level.set_attr("long_name", "sea level elevation");

  set<string> vars_to_write;
  map<string, NCSpatialVariable> list;
  surface->add_vars_to_output("big", list);
  ocean->add_vars_to_output("big", list);

  map<string, NCSpatialVariable>::iterator j = list.begin();
  while (j != list.end()) {
    vars_to_write.insert(j->first);
    ++j;
  }

  // write the states
  unsigned int record_index = 0;

  while (record_index < times.size() && times[record_index] <= grid.time->current())
    record_index++;

  ierr = nc.open(filename, PISM_WRITE, true); CHKERRQ(ierr); // append=true
  ierr = nc.def_time(grid.config.get_string("time_dimension_name"),
                     grid.config.get_string("calendar"),
                     grid.time->units()); CHKERRQ(ierr);

  while (record_index < times.size() && grid.time->current() < grid.time->end()) {

    double current_time = grid.time->current(),
      next_time = times[record_index],
      dt = next_time - current_time;

    ierr = nc.append_time(grid.config.get_string("time_dimension_name"),
                          current_time); CHKERRQ(ierr);

    char timestr[TEMPORARY_STRING_LENGTH];
    snprintf(timestr, sizeof(timestr),
             "  boundary models updated for [%s, %s] ...", 
             grid.time->date().c_str(),
             grid.time->date(next_time).c_str());
    ierr = verbPrintf(2,com,"."); CHKERRQ(ierr);
    ierr = verbPrintf(3,com,"\n%s writing result to %s ..",timestr,filename.c_str()); CHKERRQ(ierr);
    strncat(timestr,"\n",1);

    ierr = nc.append_history(timestr); CHKERRQ(ierr); // append the history

    ierr = usurf->write(nc, PISM_FLOAT); CHKERRQ(ierr);

    // update surface and ocean models' outputs:
    ierr = surface->update(current_time, dt); CHKERRQ(ierr);
    ierr = ocean->update(current_time, dt); CHKERRQ(ierr);

    ierr = surface->ice_surface_mass_flux(*climatic_mass_balance); CHKERRQ(ierr);
    ierr = surface->ice_surface_temperature(*ice_surface_temp); CHKERRQ(ierr);

    PetscReal current_sea_level;
    ierr = ocean->sea_level_elevation(current_sea_level); CHKERRQ(ierr);

    ierr = ocean->shelf_base_temperature(*shelfbasetemp); CHKERRQ(ierr);
    ierr = ocean->shelf_base_mass_flux(*shelfbasemassflux); CHKERRQ(ierr);

    sea_level.append(current_sea_level, current_time, next_time);
    sea_level.interp(current_time, next_time);

    // ask ocean and surface models to write variables:
    ierr = surface->write_variables(vars_to_write, nc); CHKERRQ(ierr);
    ierr = ocean->write_variables(vars_to_write, nc); CHKERRQ(ierr);

    // This ensures that even if a surface model wrote ice_surface_temp and climatic_mass_balance we
    // over-write them with values that were actually used by IceModel.
    ierr = climatic_mass_balance->write(nc, PISM_FLOAT); CHKERRQ(ierr);
    ierr = ice_surface_temp->write(nc, PISM_FLOAT); CHKERRQ(ierr);

    // This ensures that even if a ocean model wrote shelfbasetemp and
    // shelfbasemassflux we over-write them with values that were actually used
    // by IceModel.
    ierr = shelfbasetemp->write(nc, PISM_FLOAT); CHKERRQ(ierr);
    ierr = shelfbasemassflux->write(nc, PISM_FLOAT); CHKERRQ(ierr);

    record_index++;
    grid.time->step(dt);
  }
  ierr = verbPrintf(2,com,"\n"); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

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
    required.push_back("-times");
    ierr = show_usage_check_req_opts(com, "pclimate", required,
      "  pclimate -i IN.nc -o OUT.nc -times A:dt:B [-atmosphere <name> -surface <name> -ocean <name>] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "and set up the models:\n"
      "  -atmosphere    Chooses an atmosphere model; see User's Manual\n"
      "  -surface       Chooses a surface model; see User's Manual\n"
      "  -ocean         Chooses an ocean model; see User's Manual\n"
      ); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides, true); CHKERRQ(ierr);
    config.set("run_length_years", 0);

    bool override_used;
    ierr = PISMOptionsIsSet("-config_override", override_used); CHKERRQ(ierr);

    // set an un-documented (!) flag to limit time-steps to 1 year.
    config.set_flag("pdd_limit_timestep", true);

    IceGrid grid(com, rank, size, config);

    bool flag, times_set;
    string tmp;
    vector<double> times;
    ierr = PetscOptionsBegin(grid.com, "", "PCLIMATE options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-i", "Input file name",  inname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-times", "Specifies times to save at",
                               tmp, times_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    ierr = parse_times(grid.com, config, tmp,
                       grid.time->start(),
                       grid.time->end(),
                       times);
    if (ierr != 0) {
      PetscPrintf(grid.com, "PISM ERROR: parsing the -times argument failed.\n");
      PISMEnd();
    }

    // initialize the computational grid:
    ierr = verbPrintf(2,com, 
		      "  initializing grid from NetCDF file %s...\n", inname.c_str()); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);

    grid.time->set_start(times.front());
    grid.time->set(times.front());
    grid.time->set_end(times.back());

    mapping.init("mapping", com, rank);

    // allocate IceModelVecs needed by boundary models and put them in a dictionary:
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // read data from a PISM input file (including the projection parameters)
    PIO nc(grid, "guess_format");
    unsigned int last_record;
    bool mapping_exists;

    ierr = nc.open(inname, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_var("mapping", mapping_exists); CHKERRQ(ierr);
    ierr = nc.inq_nrecords(last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    if (mapping_exists) {
      ierr = mapping.read(inname); CHKERRQ(ierr);
      ierr = mapping.print(); CHKERRQ(ierr);
    }
    last_record -= 1;

    ierr = verbPrintf(2,com,
                      "  reading fields lat,lon,mask,thk,topg,usurf from NetCDF file %s\n"
                      "    to fill fields in PISMVars ...\n",
		      inname.c_str()); CHKERRQ(ierr);

    ierr = readIceInfoFromFile(inname, last_record, variables); CHKERRQ(ierr);

    // Initialize boundary models:
    PAFactory pa(grid, config);
    PISMAtmosphereModel *atmosphere;

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
                                outname, grid, times, mapping); CHKERRQ(ierr);

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


