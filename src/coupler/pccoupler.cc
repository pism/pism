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
  PCCDEBUG = false;  // set to true and recompile if entry and exit messages for initFromOptions(),
                     // both base class and derived classes, are needed for debugging
}


PISMClimateCoupler::~PISMClimateCoupler() { 
}


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
  //config.print(); CHKERRQ(ierr); show if desired, but IceModel already prints
  printIfDebug("ending PISMClimateCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMClimateCoupler::printIfDebug(const char *message) {
  PetscErrorCode ierr;
  if (PCCDEBUG) {  ierr = verbPrintf(1,PETSC_COMM_WORLD,message); CHKERRQ(ierr);  }
  return 0;
}


/*!
Read PISM options -i, -boot_from to determine if a PISM input or bootstrap file was given
Open the file for reading and determine its computational grid parameters; these parameters
are in returned \c LocalInterpCtx.

\e Caller of findPISMInputFile() is in charge of destroying the returned lic.
 */
PetscErrorCode PISMClimateCoupler::findPISMInputFile(char* filename, LocalInterpCtx* &lic) {
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

  // filename now contains name of PISM input file;  now check it is really there;
  // if so, read the dimensions of computational grid so that we can set up a
  // LocalInterpCtx for actual reading of climate data
  NCTool nc(grid);
  grid_info gi;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // *caller* of findPISMInputFile() is in charge of destroying
  lic = new LocalInterpCtx(gi, NULL, NULL, *grid); // 2D only

  return 0;
}


//! A virtual method which just calls specific updates.
PetscErrorCode PISMClimateCoupler::updateClimateFields(
   PetscScalar /*t_years*/, PetscScalar /*dt_years*/) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
}


//! A virtual method which writes fields associated to the derived class.
PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(PetscScalar /*t_years*/, const char */*filename*/) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
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

Derived class implementations will check user options to configure further stuff.

g->year must be valid before this can be called.
 */
PetscErrorCode PISMAtmosphereCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMAtmosphereCoupler::initFromOptions()\n");

  ierr = PISMClimateCoupler::initFromOptions(g, variables); CHKERRQ(ierr);
  
  // mean annual net ice equivalent surface mass balance rate
  ierr = vsurfmassflux.create(*g, "acab", false); CHKERRQ(ierr);
  ierr = vsurfmassflux.set_attrs(
            "",  // pism_intent is either "climate_state" or "climate_diagnostic" according
                 //    to whether this variable is kept or overwritten by a parameterization
                 //    (we don't know in the base class)
            "instantaneous net ice equivalent accumulation (ablation) rate",
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
            "",  // pism_intent is either "climate_state" or "climate_diagnostic" according
                 //    to whether this variable is kept or overwritten by a parameterization
                 //    (we don't know in the base class)
            "ice temperature at ice surface but below firn processes",
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
    ierr = dTforcing->set_dimension_units("years", ""); CHKERRQ(ierr);

    TsOffset = 0.0;
    ierr = verbPrintf(2, grid->com, 
         "  reading delta T data from forcing file %s for PISMAtmosphereCoupler...\n", dTfile); 
    CHKERRQ(ierr);
	 
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
not DOUBLE because these are expected to be imprecise at that level and not be
essential for restart accuracy.
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


//! Calls updateSurfMassFluxAndProvide() and updateSurfTempAndProvide(); ignores returned pointers.
PetscErrorCode PISMAtmosphereCoupler::updateClimateFields(
                 PetscScalar t_years, PetscScalar dt_years) {
  PetscErrorCode ierr;
  IceModelVec2* ignored;
  ierr = updateSurfMassFluxAndProvide(t_years, dt_years, ignored); CHKERRQ(ierr);
  ierr = updateSurfTempAndProvide(t_years, dt_years, ignored); CHKERRQ(ierr);
  return 0;
}


/*******************  ATMOSPHERE:  PISMConstAtmosCoupler ********************/

PISMConstAtmosCoupler::PISMConstAtmosCoupler() : PISMAtmosphereCoupler() {
  initializeFromFile = true; // default
}


//! Initializes surface mass flux and surface temperature from the PISM input file.
/*!
Because the PISMAtmosphereCoupler update procedures are not redefined, the climate
is read from file when PISM is started, but then does not change.

Non-default case: if initializeFromFile==false then nothing happens, other than
what happens in PISMAtmosphereCoupler::initFromOptions().  If you want this,
set initializeFromFile=false before calling.
 */
PetscErrorCode PISMConstAtmosCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMConstAtmosCoupler::initFromOptions()\n");

  ierr = PISMAtmosphereCoupler::initFromOptions(g, variables); CHKERRQ(ierr);

  // these values will be written into output file unchanged; they should be read
  //   as the state of the climate by future runs using the same coupler
  ierr = vsurfmassflux.set_attr("pism_intent","climate_state"); CHKERRQ(ierr);
  ierr = vsurftemp.set_attr("pism_intent","climate_state"); CHKERRQ(ierr);
  
  if (initializeFromFile) {
    char filename[PETSC_MAX_PATH_LEN];
    LocalInterpCtx* lic;
    ierr = findPISMInputFile((char*) filename, lic); CHKERRQ(ierr); // allocates lic
    ierr = verbPrintf(2, g->com,
       "  initializing constant atmospheric climate coupler ...\n"); CHKERRQ(ierr);
    ierr = verbPrintf(2, g->com, 
       "    reading net surface mass flux 'acab' from %s ...\n", filename); CHKERRQ(ierr); 
    ierr = vsurfmassflux.regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = verbPrintf(2, g->com, 
       "    reading ice surface temperature 'artm' from %s ...\n", filename); CHKERRQ(ierr); 
    ierr = vsurftemp.regrid(filename, *lic, true); CHKERRQ(ierr);
    delete lic;
  }

  printIfDebug("ending PISMConstAtmosCoupler::initFromOptions()\n");
  return 0;
}



/*******************  ATMOSPHERE:  PISMSnowModelAtmosCoupler ********************/


/* MINIMAL TEST FOR THIS SCHEME:
cd examples/eisgreen/
ncrename -v snowaccum,snowprecip green20km_y1.nc foo.nc
pcctestsm -i foo.nc -sma -pdd_std_dev 1.0 -ys 0.0 -ye 1.0 -dt 0.01 -o smamovie01.nc
pcctestsm -i foo.nc -sma -pdd_rand -pdd_std_dev 1.0 -ys 0.0 -ye 1.0 -dt 0.01 -o smarandmovie01.nc
COMPARE FIELDS acab
*/

PISMSnowModelAtmosCoupler::PISMSnowModelAtmosCoupler() : PISMAtmosphereCoupler() {
  monthlysnowtemps = NULL;
  mbscheme = NULL;
}


PISMSnowModelAtmosCoupler::~PISMSnowModelAtmosCoupler() {
  vsnowprecip.destroy();
  vsnowtemp_ma.destroy();
  vsnowtemp_mj.destroy();
  if (monthlysnowtemps != NULL) {
    delete monthlysnowtemps;
    monthlysnowtemps = NULL;
  }
  if (mbscheme != NULL) {
    delete mbscheme;
    mbscheme = NULL;
  }
}


PetscErrorCode PISMSnowModelAtmosCoupler::setLMBScheme(LocalMassBalance *usethisscheme) {
  if (mbscheme != NULL) {
    delete mbscheme;
  }
  mbscheme = usethisscheme;
  return 0;
}


PetscErrorCode PISMSnowModelAtmosCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  PetscTruth   optSet;
  printIfDebug("entering PISMSnowModelAtmosCoupler::initFromOptions()\n");

  ierr = PISMAtmosphereCoupler::initFromOptions(g, variables); CHKERRQ(ierr);

  char  filename[PETSC_MAX_PATH_LEN];
  LocalInterpCtx* lic;
  ierr = findPISMInputFile((char*)filename, lic); CHKERRQ(ierr); // allocates and initializes lic

  ierr = verbPrintf(2, g->com, 
       "  initializing atmospheric climate coupler with a snow process model ...\n"); CHKERRQ(ierr); 

  surfelev = dynamic_cast<IceModelVec2*>(variables.get("surface_altitude"));
  if (!surfelev) SETERRQ(1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2*>(variables.get("latitude"));
  if (!lat) SETERRQ(1, "ERROR: latitude is not available");

  lon = dynamic_cast<IceModelVec2*>(variables.get("longitude"));
  if (!lon) SETERRQ(1, "ERROR: longitude is not available");

  // clear out; will be overwritten by mass balance model
  ierr = vsurfmassflux.set_attr("pism_intent","climate_diagnostic"); CHKERRQ(ierr);
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr);

  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = vsnowprecip.create(*g, "snowprecip", false); CHKERRQ(ierr);
  ierr = vsnowprecip.set_attrs(
            "climate_state", 
            "mean annual ice-equivalent snow precipitation rate",
	    "m s-1", 
	    "");  // no CF standard_name ??
	    CHKERRQ(ierr);
  ierr = vsnowprecip.set_glaciological_units("m year-1");
  vsnowprecip.write_in_glaciological_units = true;

  // read snow precipitation rate from file
  ierr = verbPrintf(2, g->com, 
      "    reading mean annual ice-equivalent snow precipitation rate 'snowprecip'\n"
      "      from %s ... \n",
      filename); CHKERRQ(ierr); 
  ierr = vsnowprecip.regrid(filename, *lic, true); CHKERRQ(ierr); // fails if not found!

  // check on whether we should read monthly snow-surface temperatures from file  
  char monthlyTempsFile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-pdd_monthly_temps", 
             monthlyTempsFile, PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet == PETSC_TRUE) {
    ierr = verbPrintf(2,grid->com,
       "    reading monthly snow-surface temperatures from file %s ...\n",
       monthlyTempsFile); CHKERRQ(ierr);
    monthlysnowtemps = new MonthlyDataMaps;
    // puts month-by-month "reading ..." message to stdout:
    ierr = monthlysnowtemps->initFromFile(
             "monsnowtemp", // expects short names monsnowtemp1,...,monsnowtemp12
             "",            // choose no CF standard_name
             monthlyTempsFile, grid); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid->com,
       "    using default snow-surface temperature parameterization\n"); CHKERRQ(ierr);
    if (monthlysnowtemps != NULL)   delete monthlysnowtemps;
    monthlysnowtemps = NULL;  // test for NULL to see if using monthly temps versus parameterization
  }

  // check if user wants default or random PDD; initialize mbscheme to one of these PDDs
  if (mbscheme == NULL) { // only read user options if scheme is not chosen yet;
                          //   derived class could choose
    PetscTruth  pddRandSet, pddRepeatableSet;
    ierr = check_option("-pdd_rand", pddRandSet); CHKERRQ(ierr);
    ierr = check_option("-pdd_rand_repeatable", pddRepeatableSet); CHKERRQ(ierr);
    if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
      ierr = verbPrintf(2,grid->com,
         "    mass balance scheme chosen: PDD with simulated random daily variability ...\n");
         CHKERRQ(ierr);
      mbscheme = new PDDrandMassBalance(&config,(pddRepeatableSet == PETSC_TRUE));
    } else { // default case
      ierr = verbPrintf(2,grid->com,
         "    mass balance scheme chosen: PDD with expected value for daily variability ...\n");
         CHKERRQ(ierr);
      mbscheme = new PDDMassBalance(&config);
    }
  }
  ierr = mbscheme->init(); CHKERRQ(ierr);

  // snow temperatures used within mass balance model and from the temperature parameterization;
  //   these are ignored if monthly temps are used
  ierr = vsnowtemp_ma.create(*g, "snowtemp_ma", false); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.set_attrs(
            "climate_state",
            "mean annual snow-surface temperature used in mass balance scheme",
            "K",
            ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = vsnowtemp_ma.set(273.15); CHKERRQ(ierr);  // merely a default value
  ierr = vsnowtemp_mj.create(*g, "snowtemp_mj", false); CHKERRQ(ierr);
  ierr = vsnowtemp_mj.set_attrs(
            "climate_state",
            "mean July snow-surface temperature used in mass balance scheme",
            "K",
            ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = vsnowtemp_mj.set(273.15); CHKERRQ(ierr);  // merely a default value

  // read surface boundary condition for ice conservation of energy model from file
  //   IF AVAILABLE; derived classes may choose to reimplement
  //   PISMAtmosphereCoupler::updateSurfTempAndProvide(), which merely provides access;
  //   if read here fails, and if derived class does not reimplement, then
  //   updateSurfTempAndProvide() will fail
  ierr = verbPrintf(2, g->com, 
      "    attempting to read ice surface temperature (energy conservation boundary values) 'artm'\n"
      "      from %s ...\n",
      filename); CHKERRQ(ierr); 
  ierr = vsurftemp.regrid(filename, *lic, false); CHKERRQ(ierr); // proceeds if not found

  delete lic;

  printIfDebug("ending PISMSnowModelAtmosCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMSnowModelAtmosCoupler::writeCouplingFieldsToFile(
                              PetscScalar t_years, const char *filename) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::writeCouplingFieldsToFile(t_years,filename); CHKERRQ(ierr);

  // duplicate precipitation into file
  ierr = vsnowprecip.write(filename, NC_FLOAT); CHKERRQ(ierr);

  double dT_offset = 0.0;
  if (dTforcing != NULL)
    dT_offset = (*dTforcing)(t_years);

  // save snapshot from yearly cycle or by interpolating monthlies
  IceModelVec2 vsnowtemp;  // no work vectors, so we create a new one
  ierr = vsnowtemp.create(*grid, "snowtemp", false); CHKERRQ(ierr);
  ierr = vsnowtemp.set_attrs(
            "climate_diagnostic",
            "time-dependent snow-surface temperature used in mass balance scheme",
            "K",
            ""); CHKERRQ(ierr);  // use no CF standard_name
  PetscScalar **T;
  ierr = vsnowtemp.get_array(T);  CHKERRQ(ierr);
  if (monthlysnowtemps != NULL) {
    PetscInt curr, next;
    PetscScalar lam, **currmap, **nextmap;
    const PetscScalar t = (t_years - floor(t_years)) * secpera;
    ierr = monthlysnowtemps->getIndicesFromTime(t,curr,next,lam); CHKERRQ(ierr);
    ierr = monthlysnowtemps->vdata[curr].get_array(currmap); CHKERRQ(ierr);
    ierr = monthlysnowtemps->vdata[next].get_array(nextmap); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
      for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
        T[i][j] = monthlysnowtemps->interpolateMonthlyData(i,j,currmap,nextmap,lam) + dT_offset;
      }
    }
    ierr = monthlysnowtemps->vdata[curr].end_access(); CHKERRQ(ierr);
    ierr = monthlysnowtemps->vdata[next].end_access(); CHKERRQ(ierr);
  } else {
    PetscScalar **T_ma, **T_mj;
    const PetscScalar radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
                      sperd = 8.64e4, // exact number of seconds per day
                      julydaysec = sperd * config.get("snow_temp_july_day"),
                      t = (t_years - floor(t_years)) * secpera,
                      cosfactor = cos(radpersec * (t - julydaysec));
    ierr = vsnowtemp_ma.get_array(T_ma); CHKERRQ(ierr);
    ierr = vsnowtemp_mj.get_array(T_mj); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
      for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
        T[i][j] = T_ma[i][j] + (T_mj[i][j] - T_ma[i][j]) * cosfactor + dT_offset;
      }
    }
    ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
    ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  }
  ierr = vsnowtemp.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = vsnowtemp.destroy(); CHKERRQ(ierr);

  // no duplication of input monthly maps into output; would be done by call
  //   monthlysnowtemps->write(filename);

  return 0;
}


/*!
  The default method here is the Fausto et al parameterization scheme
  appropriate to the Greenland ice sheet.  The parameterization depends linearly
  on surface elevation, latitude, and longitude.
  
  See formulas (1) and (2) and Table 3 in \ref Faustoetal2009.
  
  Default values, stored in src/pism_config.cdl, use lines 
  'Best annual fit: This study with \f$\kappa_{\text{ma}}\f$' 
  and 'Best July fit: This study with \f$\kappa_{\text{ma}}\f$'
  from Table 3.
  
  Any scheme that has the same kind of dependence on the same quantities
  can be used just by changing the "snow_temp_fausto_..." configuration parameters 
  (in src/pism_config.cdl).  Other temperature parameterization schemes should be
  re-implement this procedure.
 */
PetscErrorCode PISMSnowModelAtmosCoupler::parameterizedUpdateSnowSurfaceTemp(
		PetscScalar /*t_years*/, PetscScalar /*dt_years*/) {
  PetscErrorCode ierr;
  const PetscScalar 
    d_ma     = config.get("snow_temp_fausto_d_ma"),      // K
    gamma_ma = config.get("snow_temp_fausto_gamma_ma"),  // K m-1
    c_ma     = config.get("snow_temp_fausto_c_ma"),      // K (degN)-1
    kappa_ma = config.get("snow_temp_fausto_kappa_ma"),  // K (degW)-1
    d_mj     = config.get("snow_temp_fausto_d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = config.get("snow_temp_fausto_gamma_mj"),
    c_mj     = config.get("snow_temp_fausto_c_mj"),
    kappa_mj = config.get("snow_temp_fausto_kappa_mj");
  
  PetscScalar **lat_degN, **lon_degE, **h, **T_ma, **T_mj;
  ierr = surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = lon->get_array(lon_degE); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);

  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      T_ma[i][j] = d_ma + gamma_ma * h[i][j] + c_ma * lat_degN[i][j] + kappa_ma * (-lon_degE[i][j]);
      T_mj[i][j] = d_mj + gamma_mj * h[i][j] + c_mj * lat_degN[i][j] + kappa_mj * (-lon_degE[i][j]);
    }
  }
  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = lon->end_access(); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  return 0;
}

//! Compute the ice surface mass flux from the snow temperature (parameterized or stored as monthly maps), a stored map of snow precipication rate, and a choice of mass balance scheme (typically a PDD).
/*!
First this method recomputes the snow temperature, either from a parameterization or from interpolating
stored monthly snow temperature maps.  More precisely, if there are no monthly temperature maps then
useTemperatureParameterizationToUpdateSnowTemp() is called .  In either case, at each point on the surface
a temperature time series with short (weekly or less) time steps is generated.
When there are no monthly temperatures, there is a parameterized yearly cycle 
of temperature and additional weather-related variability according to a normally distributed 
random temperature change for each week and grid point.

At each point on the ice surface a temperature time series, for the duration of the 
time step specified in calling this routine, is used by a LocalMassBalance object
to compute the number of positive degree days.  There are two such schemes, derived 
classes of LocalMassBalance.  There is a deterministic
(expectation) default method and random (monte carlo) method.  The standard deviation of
the temperature variation change can be controlled by option <tt>-pdd_std_dev</tt>.  
The deterministic method computes only the expected number of positive degree days
for that amount of variability \ref CalovGreve05; it is chosen by option <tt>-pdd</tt>.
The random method uses pseudo-random numbers to simulate the variability and then directly 
sums the number of positive degree days.  It is chosen by either <tt>-pdd_rand</tt> or 
<tt>-pdd_rand_repeatable</tt>.

Alternatively, if option <tt>-pdd_monthly_temps</tt> provides a NetCDF file 
with the 12 monthly temperature maps, the temperature on a given day at a given
location is found by linear interpolation (in time) of the monthly temps.

The surface mass balance is computed from the stored map of snow precipitation rate
by a call to the getMassFluxFromTemperatureTimeSeries() method of the LocalMassBalance object.
 */
PetscErrorCode PISMSnowModelAtmosCoupler::updateSurfMassFluxAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsmf) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
     t_years, dt_years, pvsmf); CHKERRQ(ierr);
 
  // set up snow temperature time series
  PetscInt     Nseries;
  ierr = mbscheme->getNForTemperatureSeries(
             t_years * secpera, dt_years * secpera, Nseries); CHKERRQ(ierr);
  PetscScalar  *Tseries = new PetscScalar[Nseries];

  // times at which we compute snow temps are
  //    tseries, tseries + dtseries, ..., tseries + (Nseries-1) * dtseries;
  const PetscScalar tseries = (t_years - floor(t_years)) * secpera,
                    dtseries = (dt_years * secpera) / ((PetscScalar) Nseries);
  
  // constants related to standard yearly cycle
  const PetscScalar
     radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
     sperd = 8.64e4, // exact number of seconds per day
     julydaysec = sperd * config.get("snow_temp_july_day");
  
  // get access to snow temperatures (and update parameterization if appropriate)
  PetscScalar **T_ma, **T_mj, **snowrate, **lat_degN, **smflux;
  if (monthlysnowtemps == NULL) {
    ierr = parameterizedUpdateSnowSurfaceTemp(t_years,dt_years); CHKERRQ(ierr);
  }
  ierr = vsnowtemp_ma.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);
  ierr = vsnowprecip.get_array(snowrate);  CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = vsurfmassflux.get_array(smflux);  CHKERRQ(ierr);

  // run through grid
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {

      // build temperature time series at point i,j
      for (PetscInt k = 0; k < Nseries; ++k) {
        const PetscScalar  tk = tseries + ((PetscScalar) k) * dtseries;
        if (monthlysnowtemps != NULL) {
          // interpolate monthly temps; not a great memory access pattern ...
          PetscInt curr, next;
          PetscScalar lam;
          ierr = monthlysnowtemps->getIndicesFromTime(tk,curr,next,lam); CHKERRQ(ierr);
          PetscScalar **currmap, **nextmap;
          ierr = monthlysnowtemps->vdata[curr].get_array(currmap); CHKERRQ(ierr);
          ierr = monthlysnowtemps->vdata[next].get_array(nextmap); CHKERRQ(ierr);
          Tseries[k] = monthlysnowtemps->interpolateMonthlyData(i,j,currmap,nextmap,lam);
          ierr = monthlysnowtemps->vdata[curr].end_access(); CHKERRQ(ierr);
          ierr = monthlysnowtemps->vdata[next].end_access(); CHKERRQ(ierr);
        } else { // default case
          // use corrected formula (4) in \ref Faustoetal2009, the standard yearly cycle
          Tseries[k] = T_ma[i][j] + (T_mj[i][j] - T_ma[i][j]) * cos(radpersec * (tk - julydaysec));
        }
	
	// apply the temperature offset if needed
	if (dTforcing != NULL)
	  Tseries[k] += (*dTforcing)(t_years + k * dt_years / Nseries);

      }

      // if we can, set mass balance parameters according to formula (6) in
      //   \ref Faustoetal2009
      PDDMassBalance* pddscheme = dynamic_cast<PDDMassBalance*>(mbscheme);
      if (pddscheme != NULL) {
        ierr = pddscheme->setDegreeDayFactorsFromSpecialInfo(lat_degN[i][j],T_mj[i][j]); CHKERRQ(ierr);
      }

      
      // get surface mass balance at point i,j
      smflux[i][j] = mbscheme->getMassFluxFromTemperatureTimeSeries(
                        tseries, dtseries, Tseries, Nseries,snowrate[i][j]);  // units of day^-1 K-1
    }
  }
  
  ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  ierr = vsnowprecip.end_access();  CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = vsurfmassflux.end_access(); CHKERRQ(ierr);

  delete [] Tseries;
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

