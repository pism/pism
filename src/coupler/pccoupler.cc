// Copyright (C) 2008-2009 Ed Bueler, Ricarda Winkelmann and Constantine Khroulev
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


#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModelVec.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/nc_util.hh"
#include "../base/NCVariable.hh"
#include "../base/Timeseries.hh"
#include "pccoupler.hh"
// we do NOT depend on IceModel.hh; this is deliberate!


/******************* VIRTUAL BASE CLASS:  PISMClimateCoupler ********************/

PISMClimateCoupler::PISMClimateCoupler() {
  grid = NULL;
  PCCDEBUG = false;  // set this to true and recompile if entry and exit messages
                     //   for initFromOptions() are needed for debugging
}


PISMClimateCoupler::~PISMClimateCoupler() { 
}


//! Caller of this routine provides an IceGrid to the coupler.
/*!
The IceGrid *g is used to set up access to the config file.
 */
PetscErrorCode PISMClimateCoupler::initFromOptions(IceGrid* g, const PISMVars &/*variables*/) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMClimateCoupler::initFromOptions()\n");

  grid = g;
  config.init("pism_config", grid->com, grid->rank);
  char alt_config[PETSC_MAX_PATH_LEN];
  PetscTruth use_alt_config;
  ierr = PetscOptionsGetString(PETSC_NULL, "-config", alt_config, 
                               PETSC_MAX_PATH_LEN, &use_alt_config); CHKERRQ(ierr);
  if (use_alt_config) {
    ierr = config.read(alt_config); CHKERRQ(ierr);
  } else {
    ierr = config.read(PISM_DefaultConfigFile); CHKERRQ(ierr);
  }

  printIfDebug("ending PISMClimateCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMClimateCoupler::printIfDebug(const char *message) {
  PetscErrorCode ierr;
  if (PCCDEBUG) {  ierr = verbPrintf(1,PETSC_COMM_WORLD,message); CHKERRQ(ierr);  }
  return 0;
}


  // 


//! Provides access, from the coupler, to the main PISM input file.
/*!
Since climate fields may be in the same file as the one used by IceModel for input,
we need to find it and get info needed for reading it.  Thus this method is normally
a helper routine used by derived classes of PISMClimateCoupler.

Reads PISM options -i, -boot_from to determine if a PISM input or bootstrap file was given.
Opens the file for reading and determine its computational grid parameters.  These parameters
are in returned \c LocalInterpCtx.

\e Caller of findPISMInputFile() is in charge of destroying the returned lic.

Sets lic to NULL if the input file was a -i input file.  This case means that
we don't need to regrid, and we can just read data straight.  In this same case the
returned flag \c regrid is false, but the returned integer \c start gives the 
time dimension index to use in reading.

If -boot_from was used, returns the \c LocalInterpCtx needed to read it, and
\c regrid is true.  In this case the returned integer \c start should be ignored. 

Argument filename needs to be pre-allocated.
 */
PetscErrorCode PISMClimateCoupler::findPISMInputFile(char* filename, LocalInterpCtx* &lic,
						     bool &regrid, int &start) {
  PetscErrorCode ierr;
  PetscTruth i_set, boot_from_set;

  if (grid == NULL) {  SETERRQ(1,"findPISMInputFile(): grid not initialized");  }

  // read file names:
  char i_file[PETSC_MAX_PATH_LEN], boot_from_file[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-i", i_file, 
			       PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from", boot_from_file, 
			       PETSC_MAX_PATH_LEN, &boot_from_set); CHKERRQ(ierr);
  if (i_set) {
    if (boot_from_set) {
      ierr = PetscPrintf(grid->com,
	"PISMClimateCoupler ERROR: both '-i' and '-boot_from' are used. Exiting...\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    strcpy(filename, i_file);
  }
  else if (boot_from_set) {
    strcpy(filename, boot_from_file);
  }

  // filename now contains name of PISM input (bootstrapping) file; now check
  // it is really there; if so, read the dimensions of computational grid so
  // that we can set up a LocalInterpCtx for actual reading of climate data
  NCTool nc(grid);
  int last_record;
  grid_info gi;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
  last_record -= 1;
  ierr = nc.close(); CHKERRQ(ierr);

  if (boot_from_set) {
    // *caller* of findPISMInputFile() is in charge of destroying
    lic = new LocalInterpCtx(gi, NULL, NULL, *grid); // 2D only
    regrid = true;
    start = 0;
  } else {
    lic = NULL;
    regrid = false;
    start = last_record;
  }

  return 0;
}


PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(
            PetscScalar /*t_years*/, const char */*filename*/) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and is not implemented");
}


//! A virtual method which just calls specific updates.
PetscErrorCode PISMClimateCoupler::updateClimateFields(
            PetscScalar /*t_years*/, PetscScalar /*dt_years*/) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and is not implemented");
}

//! Sets dt_years to min(dt_years, maximum step the coupler can take).
/*!
  Does nothing in the base class.
 */
PetscErrorCode PISMClimateCoupler::max_timestep(PetscScalar /*t_years*/, PetscScalar &/*dt_years*/) {
  return 0;
}

/******************* ATMOSPHERE:  PISMAtmosphereCoupler ********************/

PISMAtmosphereCoupler::PISMAtmosphereCoupler() : PISMClimateCoupler() {
  dTforcing = PETSC_NULL;
  TsOffset = 0.0;
}


PISMAtmosphereCoupler::~PISMAtmosphereCoupler() {
  vsurfmassflux.destroy();  // destroy if created
  vsurftemp.destroy();
  if (dTforcing != PETSC_NULL) {
    delete dTforcing; // calls destructor for this Timeseries instance
    dTforcing = PETSC_NULL;
  }
  TsOffset = 0.0;
}


//! Initialize a PISMAtmosphereCoupler by allocating space for surface mass flux and surface temperature variables.
/*!
Allocates space and sets attributes, including CF standard_name, for the two essential fields,
namely the two fields to which IceModel needs access.

The short names "acab" and "artm" for these two fields match GLIMMER (& CISM, presumably).

Derived class implementations of this routine may check user options to configure
further stuff.

g->year must be valid before this can be called.  (FIXME:  why?)

If option -dTforcing is used, we find and read the temperature offsets from the given
file into memory, by calling Timeseries.read() for Timeseries* dTforcing.
 */
PetscErrorCode PISMAtmosphereCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMAtmosphereCoupler::initFromOptions()\n");

  ierr = PISMClimateCoupler::initFromOptions(g, variables); CHKERRQ(ierr);
  
  // mean annual net ice equivalent surface mass balance rate
  ierr = vsurfmassflux.create(*g, "acab", false); CHKERRQ(ierr);
  ierr = vsurfmassflux.set_attrs(
            "",  // pism_intent is either "climate_state" or "climate_diagnostic"
                 //    according to whether this variable is kept or overwritten
                 //    by a parameterization (we don't know in the base class)
            "instantaneous ice-equivalent accumulation (ablation) rate",
	    "m s-1",  // m *ice-equivalent* per second
	    "land_ice_surface_specific_mass_balance");  // CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfmassflux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vsurfmassflux.write_in_glaciological_units = true;
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr);  // merely a default value

  // annual mean air temperature at "ice surface", at level below all firn processes
  //   (e.g. "10 m" ice temperatures)
  ierr = vsurftemp.create(*g, "artm", false); CHKERRQ(ierr);
  ierr = vsurftemp.set_attrs(
            "",  // pism_intent is either "climate_state" or "climate_diagnostic"
                 //    according to whether this variable is kept or overwritten
                 //    by a parameterization (we don't know in the base class)
            "annual average ice temperature at ice surface but below firn processes",
            "K", 
            "");  // PROPOSED CF standard_name = land_ice_surface_temperature_below_firn
            CHKERRQ(ierr);
  ierr = vsurftemp.set(273.15); CHKERRQ(ierr);  // merely a default value

  // check user option -dTforcing for a surface temperature forcing data set
  char dTfile[PETSC_MAX_PATH_LEN];
  PetscTruth dTforceSet;
  if (dTforcing != PETSC_NULL) {
    SETERRQ(1, "dTforcing!=PETSC_NULL in PISMAtmosphereCoupler::initFromOptions()\n");
  }
  ierr = PetscOptionsGetString(PETSC_NULL, "-dTforcing", dTfile,
                               PETSC_MAX_PATH_LEN, &dTforceSet); CHKERRQ(ierr);
  if (dTforceSet == PETSC_TRUE) {
    dTforcing = new Timeseries(grid, "delta_T", "t");
    ierr = dTforcing->set_units("Celsius", ""); CHKERRQ(ierr);
    ierr = dTforcing->set_attr("long_name", "surface temperature offsets"); CHKERRQ(ierr);
    ierr = dTforcing->set_dimension_units("years", ""); CHKERRQ(ierr);

    TsOffset = 0.0;
    ierr = verbPrintf(2, grid->com, 
         "  reading delta T data from forcing file %s for PISMAtmosphereCoupler...\n",
         dTfile);  CHKERRQ(ierr);
	 
    ierr = dTforcing->read(dTfile); CHKERRQ(ierr);
  }

  printIfDebug("ending PISMAtmosphereCoupler::initFromOptions()\n");
  return 0;
}


//! Writes surface mass flux and surface temperature to prepared file.
/*!
Assumes file is prepared in the sense that it exists and that global attributes are
already written.  See IceModel::dumpToFile() for how main PISM output file is
prepared.  Calls here do handle opening and closing the file.  We write in FLOAT 
not DOUBLE because these are expected to be imprecise at that level and they will not be
used for restarting.
 */
PetscErrorCode PISMAtmosphereCoupler::writeCouplingFieldsToFile(
		PetscScalar /*t_years*/, const char *filename) {
  PetscErrorCode ierr;
  
  ierr = vsurfmassflux.write(filename, NC_FLOAT); CHKERRQ(ierr);
  if (dTforcing != PETSC_NULL) {
    ierr = vsurftemp.shift(-TsOffset); CHKERRQ(ierr); // return to unshifted state
    ierr = vsurftemp.write(filename, NC_FLOAT); CHKERRQ(ierr);
    ierr = vsurftemp.shift(TsOffset); CHKERRQ(ierr);  // re-apply the offset
  } else {
    ierr = vsurftemp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  }

  // also append to surface temperature offset time series
  NCTool nc(grid);
  bool variable_exists;
  ierr = nc.open_for_writing(filename, true, true);
  // append == true, check_dims == true
  ierr = nc.find_variable("surftempoffset", NULL, variable_exists); CHKERRQ(ierr);
  if (!variable_exists) {
    ierr = nc.create_timeseries("surftempoffset", "surface temperature offset",
				"Celsius", NC_FLOAT, NULL);
    CHKERRQ(ierr);
  }
  ierr = nc.append_timeseries("surftempoffset", TsOffset); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}


//! Calls updateSurfMassFluxAndProvide() and updateSurfTempAndProvide(); ignors returned pointers.
PetscErrorCode PISMAtmosphereCoupler::updateClimateFields(
                 PetscScalar t_years, PetscScalar dt_years) {
  PetscErrorCode ierr;
  IceModelVec2* ignored;
  ierr = updateSurfMassFluxAndProvide(t_years, dt_years, ignored); CHKERRQ(ierr);
  ierr = updateSurfTempAndProvide(t_years, dt_years, ignored); CHKERRQ(ierr);
  return 0;
}


//! Provides access to vsurfmassflux.  No update of vsurfmassflux.  Derived class versions generally will update.
PetscErrorCode PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
    PetscScalar /*t_years*/, PetscScalar /*dt_years*/,
    IceModelVec2* &pvsmf) {
  if (vsurfmassflux.was_created())
    pvsmf = &vsurfmassflux;
  else {  SETERRQ(1,"vsurfmassflux not created in updateSurfMassFluxAndProvide()");  }
  return 0;
}


//! Updates vsurftemp using -dTforcing (if it is on) and provides access to vsurftemp.  Derived class versions may do more updating.
PetscErrorCode PISMAtmosphereCoupler::updateSurfTempAndProvide(
    PetscScalar t_years, PetscScalar /*dt_years*/,
    IceModelVec2* &pvst) {
  PetscErrorCode ierr;
  if (vsurftemp.was_created())
    pvst = &vsurftemp;
  else {  SETERRQ(1,"vsurftemp not created in updateSurfTempAndProvide()");  }

  if (dTforcing != PETSC_NULL) {
    double old_offset = TsOffset;
    TsOffset = (*dTforcing)(t_years);

    ierr = verbPrintf(5,grid->com,
       "PISMAtmosphereCoupler says: read TsOffset=%.6f from -dTforcing data\n",
       TsOffset); CHKERRQ(ierr);
    ierr = vsurftemp.shift(TsOffset - old_offset); CHKERRQ(ierr);  // apply the new offset
  }

  return 0;
}


/*******************  ATMOSPHERE:  PISMConstAtmosCoupler ********************/

PISMConstAtmosCoupler::PISMConstAtmosCoupler() : PISMAtmosphereCoupler() {
  initializeFromFile = true; // default
}


//! Normally initializes surface mass flux and surface temperature from the PISM input file.
/*!
Default case: initializeFromFile==true then reads vsurfmassflux and vsurftemp
from file.

Non-default case: if initializeFromFile==false then nothing happens, other than
what happens in PISMAtmosphereCoupler::initFromOptions().  If you want this,
set initializeFromFile=false before calling.

In either case, because the PISMAtmosphereCoupler update procedures are not
redefined, the climate does not change and is the same is it is when initialized. 
 */
PetscErrorCode PISMConstAtmosCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMConstAtmosCoupler::initFromOptions()\n");

  ierr = PISMAtmosphereCoupler::initFromOptions(g, variables); CHKERRQ(ierr);

  // because these values will be written into output file unchanged,
  //    they are the state of the climate for future runs using the same coupler
  ierr = vsurfmassflux.set_attr("pism_intent","climate_state"); CHKERRQ(ierr);
  ierr = vsurftemp.set_attr("pism_intent","climate_state"); CHKERRQ(ierr);
  
  if (initializeFromFile) {
    char filename[PETSC_MAX_PATH_LEN];
    LocalInterpCtx* lic;
    bool regrid;
    int start;
    // next call allocates lic (if necessary):
    ierr = findPISMInputFile((char*) filename, lic, regrid, start); CHKERRQ(ierr);
    ierr = verbPrintf(2, g->com,
       "  initializing constant atmospheric climate coupler ...\n"); CHKERRQ(ierr);

    ierr = verbPrintf(2, g->com, 
       "    reading net surface mass flux 'acab' from %s ...\n", filename); CHKERRQ(ierr); 
    if (regrid) {
      ierr = vsurfmassflux.regrid(filename, *lic, true); CHKERRQ(ierr);
    } else {
      ierr = vsurfmassflux.read(filename, start); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2, g->com, 
       "    reading ice surface temperature 'artm' from %s ...\n", filename); CHKERRQ(ierr); 
    if (regrid) {
      ierr = vsurftemp.regrid(filename, *lic, true); CHKERRQ(ierr);
    } else {
      ierr = vsurftemp.read(filename, start); CHKERRQ(ierr);
    }

    delete lic;			// deleting a NULL pointer is OK
  }

  printIfDebug("ending PISMConstAtmosCoupler::initFromOptions()\n");
  return 0;
}


/*******************  OCEAN:  PISMOceanCoupler ********************/

PISMOceanCoupler::PISMOceanCoupler() : PISMClimateCoupler() {
  reportInitializationToStdOut = true;  // derived classes etc. can turn off before calling
                                        // initFromOptions(), but it's on by default
  dSLforcing = PETSC_NULL;
  seaLevel = 0.0; // the obvious default value
}


PISMOceanCoupler::~PISMOceanCoupler() {
  vshelfbasetemp.destroy();
  vshelfbasemassflux.destroy();
  if (dSLforcing != PETSC_NULL) {
    seaLevel = 0.0;
    delete dSLforcing;
    dSLforcing = PETSC_NULL;
  }
}


/*!
Derived class implementations will check user options to configure the PISMOceanCoupler.
This version allocates space and sets attributes for the two essential fields.

g->year must be valid before this can be called.
 */
PetscErrorCode PISMOceanCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMOceanCoupler::initFromOptions()\n");

  ierr = PISMClimateCoupler::initFromOptions(g, variables); CHKERRQ(ierr);

  // no report to std out; otherwise reportInitializationToStdOut should be checked before report

  // ice boundary tempature at the base of the ice shelf
  ierr = vshelfbasetemp.create(*g, "shelfbtemp", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = vshelfbasetemp.set_attrs(
           "climate_state", "absolute temperature at ice shelf base",
	   "K", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_temperature
  ierr = vshelfbasetemp.set(273.15); CHKERRQ(ierr);  // merely a default value to clear nonsense

  // ice mass balance rate at the base of the ice shelf; sign convention for vshelfbasemass
  //   matches standard sign convention for basal melt rate of grounded ice
  ierr = vshelfbasemassflux.create(*g, "shelfbmassflux", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = vshelfbasemassflux.set_attrs(
           "climate_state", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
	   "m s-1", ""); CHKERRQ(ierr); 
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  ierr = vshelfbasemassflux.set(0.0); CHKERRQ(ierr);  // merely a default value to clear nonsense
  // rescales from m/s to m/a when writing to NetCDF and std out:
  vshelfbasemassflux.write_in_glaciological_units = true;
  ierr = vshelfbasemassflux.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  char dSLfile[PETSC_MAX_PATH_LEN];
  PetscTruth dSLforceSet;
  if (dSLforcing != PETSC_NULL) {
    SETERRQ(2, "dSLforcing should be PETSC_NULL at start of PISMOceanCoupler::initFromOptions()\n");
  }
  ierr = PetscOptionsGetString(PETSC_NULL, "-dSLforcing", dSLfile,
                               PETSC_MAX_PATH_LEN, &dSLforceSet); CHKERRQ(ierr);

  if (dSLforceSet == PETSC_TRUE) {
    dSLforcing = new Timeseries(grid, "delta_sea_level", "t");
    ierr = dSLforcing->set_units("m", ""); CHKERRQ(ierr);
    ierr = dSLforcing->set_dimension_units("years", ""); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid->com, 
		      "  reading delta sea level data from forcing file %s ...\n", 
		      dSLfile); CHKERRQ(ierr);
    ierr = dSLforcing->read(dSLfile); CHKERRQ(ierr);
  }

  printIfDebug("ending PISMOceanCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMOceanCoupler::writeCouplingFieldsToFile(
		PetscScalar /*t_years*/, const char *filename) {
  PetscErrorCode ierr;
  
  // We assume file is prepared in the sense that it exists and that global attributes 
  //   are already written.  See IceModel::dumpToFile() for how main PISM output file is
  //   prepared.  Note calls here handle opening and closing the file.  We write in
  //   FLOAT not DOUBLE because these are expected to be for diagnosis, not restart etc.
  ierr = vshelfbasetemp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = vshelfbasemassflux.write(filename, NC_FLOAT); CHKERRQ(ierr);

  NCTool nc(grid);
  bool variable_exists;
  ierr = nc.open_for_writing(filename, true, true);
  // append == true, check_dims == true
  ierr = nc.find_variable("sealevel", NULL, variable_exists); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = nc.create_timeseries("sealevel", "sea level", "meters", NC_FLOAT, NULL);
    CHKERRQ(ierr);
  }

  ierr = nc.append_timeseries("sealevel", seaLevel); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}


//! Updates all the ocean climate fields.
PetscErrorCode PISMOceanCoupler::updateClimateFields(
        PetscScalar t_years, PetscScalar dt_years) {
  PetscErrorCode ierr;
  IceModelVec2* ignored;
  ierr = updateShelfBaseMassFluxAndProvide(t_years, dt_years, ignored); CHKERRQ(ierr);
  ierr = updateShelfBaseTempAndProvide(t_years, dt_years, ignored); CHKERRQ(ierr);
  ierr = updateSeaLevelElevation(t_years, dt_years, NULL); CHKERRQ(ierr);
  return 0;
}


//! Provides access to vshelfbasemassflux.  No update of vshelfbasemassflux.  Real ocean models in derived classes will update.
PetscErrorCode PISMOceanCoupler::updateShelfBaseMassFluxAndProvide(
       PetscScalar /*t_years*/, PetscScalar /*dt_years*/,
       IceModelVec2* &pvsbmf) {

  if (vshelfbasemassflux.was_created())
    pvsbmf = &vshelfbasemassflux;
  else {  
    SETERRQ(1,"vvshelfbasemassflux not created in PISMOceanCoupler::updatehelfBaseMassFluxAndProvide()");
  }

  return 0;
}


//! Provides access to vshelfbasetemp.  No update of vshelfbasetemp.  Real ocean models in derived classes will update.
PetscErrorCode PISMOceanCoupler::updateShelfBaseTempAndProvide(
         PetscScalar /*t_years*/, PetscScalar /*dt_years*/,
         IceModelVec2* &pvsbt) {
  // printIfDebug("entering PISMOceanCoupler::updateShelfBaseTempAndProvide()\n");

  if (vshelfbasetemp.was_created())
    pvsbt = &vshelfbasetemp;
  else {  
    SETERRQ(1,"vvshelfbasetemp not created in PISMOceanCoupler::updatehelfBaseTempAndProvide()");
  }

  //  printIfDebug("leaving PISMOceanCoupler::updateShelfBaseTempAndProvide()\n");
  return 0;
}


//! Updates the sea level (using -dSLforcing, if it is on) and sets \c new_sea_level (if not NULL).
PetscErrorCode PISMOceanCoupler::updateSeaLevelElevation(PetscReal t_years, PetscReal /*dt_years*/,
							 PetscReal *new_sea_level) {
  PetscErrorCode ierr;

  if (dSLforcing != PETSC_NULL) {
    // read the new sea level (delta from modern)
    seaLevel = (*dSLforcing)(t_years);

    ierr = verbPrintf(5,grid->com,"read newSeaLevel=%.6f from -dSLforcing climate data\n",
       seaLevel); CHKERRQ(ierr);
    // comment: IceModel::updateSurfaceElevationAndMask() needs to be called
    // before effect of changed sea level is seen in ice dynamics (e.g. on
    // grounding line)
  }
  
  if (new_sea_level != NULL) *new_sea_level = seaLevel;

  return 0;
}


/*******************  OCEAN:  PISMConstOceanCoupler ********************/

/* historical note only:  PISMConstOceanCoupler is the only (quite simplified!) 
ocean model provided with the PISM source.  It is based on a 
nontrivial PISMClimateCoupler for ice shelves
which was written by Ricarda Winkelmann (PIK), with input from Bueler, on 2 Dec 2008.
That source is in icerepo/Potsdam-Antarctica, a password-protected svn repo.
Contact Bueler for further info on that. */

PISMConstOceanCoupler::PISMConstOceanCoupler() : PISMOceanCoupler() {
  constOceanHeatFlux = 0.5;   // W m-2 = J m-2 s-1; naively chosen default value
        // default value possibly irrelevant as long as it is pretty small;
        //   0.5 W m-2 is about 4 times more heating than peak of Shapiro&Ritzwoller (2004)
        //   geothermal fluxes for Antarctica of about 130 mW/m^2;
        //   0.5 W m-2 yields  0.051758 m a-1 = 5.2 cm a-1  basal melt rate, ice 
        //   thickness per time, in updateShelfBaseMassFluxAndProvide() below
        
        // alternative: a rate of zero might do no harm; note heat flux immediately
        //   becomes a basal net mass balance;
        //   Lingle et al (1991; "A flow band model of the Ross Ice Shelf ..."
        //   JGR 96 (B4), pp 6849--6871) gives 0.02 m/a freeze-on at one point as only 
        //   measurement available at that time (one ice core) and also gives
        //   0.16 m/a melting as average rate necessary to maintain equilibrium,
        //   but points out variability in -0.5 m/a (i.e. melting) to 
        //   +1.0 m/a (freeze-on) range from a flow band model (figure 5)
}


PetscErrorCode PISMConstOceanCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMConstOceanCoupler::initFromOptions()\n");

  ierr = PISMOceanCoupler::initFromOptions(g, variables); CHKERRQ(ierr);

  thk = dynamic_cast<IceModelVec2*>(variables.get("land_ice_thickness"));
  if (!thk) SETERRQ(1, "ERROR: land_ice_thickness is not available");
  
  if (reportInitializationToStdOut) {
    ierr = verbPrintf(2, g->com, 
       "  initializing constant sub-ice shelf ocean climate:\n"
       "    heat flux from ocean set to %.3f W m-2 (determines mass balance)\n"
       "    ice shelf base temperature set to pressure-melting temperature\n",
       constOceanHeatFlux); CHKERRQ(ierr); 
  }

  printIfDebug("ending PISMConstOceanCoupler::initFromOptions()\n");
  return 0;
}


//! By default, does not write vshelfbasetemp and vshelfbasemassflux fields.
PetscErrorCode PISMConstOceanCoupler::writeCouplingFieldsToFile(
		PetscScalar /*t_years*/, const char *filename) {
  PetscErrorCode ierr;
  NCTool nc(grid);
  bool variable_exists;
  ierr = nc.open_for_writing(filename, true, true);
  // append == true, check_dims == true
  ierr = nc.find_variable("sealevel", NULL, variable_exists); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = nc.create_timeseries("sealevel", "sea level", "meters", NC_FLOAT, NULL);
    CHKERRQ(ierr);
  }

  ierr = nc.append_timeseries("sealevel", seaLevel); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMConstOceanCoupler::updateShelfBaseTempAndProvide(
                  PetscScalar /*t_years*/, PetscScalar /*dt_years*/,
                  IceModelVec2* &pvsbt) {
  PetscErrorCode ierr;

  // ignores everything from IceModel except ice thickness; also ignores t_years and dt_years
  const PetscScalar icerho       = 910.0,   // kg/m^3   ice shelf mean density
                    oceanrho     = 1028.0,  // kg/m^3   sea water mean density
                    beta_CC_grad = 8.66e-4, // K/m      Clausius-Clapeyron gradient
                    T0           = 273.15;  // K        triple point = naively assumed
                                            //          sea water temperature at sea level
  PetscScalar **H, **temp;
  ierr = thk->get_array(H);   CHKERRQ(ierr);
  ierr = vshelfbasetemp.get_array (temp); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      const PetscScalar shelfbaseelev = - (icerho / oceanrho) * H[i][j];
      // temp is set to melting point at depth
      temp[i][j] = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so this is below T0
    }
  }
  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = vshelfbasetemp.end_access(); CHKERRQ(ierr);
  
  pvsbt = &vshelfbasetemp;
  return 0;                                 
}


PetscErrorCode PISMConstOceanCoupler::updateShelfBaseMassFluxAndProvide(
         PetscScalar /*t_years*/, PetscScalar /*dt_years*/,
         IceModelVec2* &pvsbmf) {
  PetscErrorCode ierr;

  const PetscScalar icelatentheat = 3.35e5,   // J kg-1   ice latent heat capacity
                    icerho        = 910.0,    // kg m-3   ice shelf mean density
                    // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
                    meltrate      = constOceanHeatFlux / (icelatentheat * icerho); // m s-1

  // vshelfbasemassflux is positive if ice is melting (flux into ocean);
  //    see metadata set in PISMOceanCoupler::initFromOptions()
  ierr = vshelfbasemassflux.set(meltrate); CHKERRQ(ierr);

  pvsbmf = &vshelfbasemassflux;
  return 0;
}

