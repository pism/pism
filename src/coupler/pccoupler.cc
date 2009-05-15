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


PetscErrorCode PISMClimateCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMClimateCoupler::initFromOptions()\n");
  grid = g;
  config.init("pism_config", *grid);
  char alt_config[PETSC_MAX_PATH_LEN];
  PetscTruth use_alt_config;
  ierr = PetscOptionsGetString(PETSC_NULL, "-config", alt_config, 
                               PETSC_MAX_PATH_LEN, &use_alt_config); CHKERRQ(ierr);
  if (use_alt_config) {
    ierr = config.read(alt_config); CHKERRQ(ierr);
  } else {
    ierr = config.read(PISM_DEFAULT_CONFIG_FILE); CHKERRQ(ierr);
  }
  config.print(); CHKERRQ(ierr); // FIXME: desired?
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
  const PetscScalar /*t_years*/, const PetscScalar /*dt_years*/, void */*iceInfoNeeded*/) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
}


//! A virtual method which writes fields associated to the derived class.
PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(const char */*filename*/) {
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
    TsOffset = 0.0;
    delete dTforcing; // calls destructor for this PISMClimateForcing instance
    dTforcing = PETSC_NULL;
  }
}


//! Initialize a PISMAtmosphereCoupler by allocating space for surface mass flux and surface temperature variables.
/*!
Allocates space and sets attributes, including CF standard_name, for the two essential fields,
namely the two fields to which IceModel needs access.

The short names "acab" and "artm" for these two fields match GLIMMER (& CISM, presumably).

Derived class implementations will check user options to configure further stuff.

g->year must be valid before this can be called.
 */
PetscErrorCode PISMAtmosphereCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMAtmosphereCoupler::initFromOptions()\n");

  ierr = PISMClimateCoupler::initFromOptions(g); CHKERRQ(ierr);
  
  // mean annual net ice equivalent surface mass balance rate
  ierr = vsurfmassflux.create(*g, "acab", false); CHKERRQ(ierr);
  ierr = vsurfmassflux.set_attrs(
            "climate_state", 
            "instantaneous net ice equivalent accumulation (ablation) rate",
	    "m s-1",  // m *ice-equivalent* per second
	    "land_ice_surface_specific_mass_balance");  // CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfmassflux.set_glaciological_units("m year-1");
  vsurfmassflux.write_in_glaciological_units = true;
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr);  // merely a default value

  // annual mean air temperature at "ice surface", at level below all firn processes
  //   (e.g. "10 m" ice temperatures)
  ierr = vsurftemp.create(*g, "artm", false); CHKERRQ(ierr);
  ierr = vsurftemp.set_attrs(
            "climate_state",
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
    dTforcing = new PISMClimateForcing;
    TsOffset = 0.0;
    int stat, ncid = 0;
    ierr = verbPrintf(2, grid->com, 
         "  reading delta T data from forcing file %s for PISMAtmosphereCoupler...\n", dTfile); 
         CHKERRQ(ierr);
    if (grid->rank == 0) {   stat = nc_open(dTfile, 0, &ncid); CHKERRQ(nc_check(stat));   }
    ierr = dTforcing->readClimateForcingData(grid->com, grid->rank, ncid, 
         grid->year,PCF_DELTA_T); CHKERRQ(ierr);
    if (grid->rank == 0) {   stat = nc_close(ncid); CHKERRQ(nc_check(stat));  }
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
PetscErrorCode PISMAtmosphereCoupler::writeCouplingFieldsToFile(const char *filename) {
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
				"degree_Celsius", NC_FLOAT, NULL);
    CHKERRQ(ierr);
  }
  ierr = nc.append_timeseries("surftempoffset", TsOffset); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}


//! Provides access to vsurfmassflux.  No update of vsurfmassflux.  Derived class versions generally will update.
PetscErrorCode PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
    const PetscScalar /*t_years*/, const PetscScalar /*dt_years*/, 
    void */*iceInfoNeeded*/, IceModelVec2* &pvsmf) {
  if (vsurfmassflux.was_created())
    pvsmf = &vsurfmassflux;
  else {  SETERRQ(1,"vsurfmassflux not created in updateSurfMassFluxAndProvide()");  }
  return 0;
}


//! Updates vsurftemp using -dTforcing (if it is on) and provides access to vsurftemp.  Derived class versions may do more updating.
PetscErrorCode PISMAtmosphereCoupler::updateSurfTempAndProvide(
  const PetscScalar t_years, const PetscScalar /*dt_years*/, 
  void */*iceInfoNeeded*/, IceModelVec2* &pvst) {
  PetscErrorCode ierr;

  if (vsurftemp.was_created())
    pvst = &vsurftemp;
  else {  SETERRQ(1,"vsurftemp not created in updateSurfTempAndProvide()");  }

  if (dTforcing != PETSC_NULL) {
    ierr = vsurftemp.shift(-TsOffset); CHKERRQ(ierr); // return to unshifted state
    ierr = dTforcing->updateFromClimateForcingData(t_years,&TsOffset); CHKERRQ(ierr); // read a new offset
    ierr = verbPrintf(5,grid->com,
       "PISMAtmosphereCoupler says: read TsOffset=%.6f from -dTforcing data\n",
       TsOffset); CHKERRQ(ierr);
    ierr = vsurftemp.shift(TsOffset); CHKERRQ(ierr);  // apply the offset
  }

  return 0;
}


//! Calls updateSurfMassFluxAndProvide() and updateSurfTempAndProvide(), but ignores returned pointers.
PetscErrorCode PISMAtmosphereCoupler::updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded) {
  PetscErrorCode ierr;
  IceModelVec2* ignored;
  ierr = updateSurfMassFluxAndProvide(t_years, dt_years, iceInfoNeeded, ignored); CHKERRQ(ierr);
  ierr = updateSurfTempAndProvide(t_years, dt_years, iceInfoNeeded, ignored); CHKERRQ(ierr);
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
PetscErrorCode PISMConstAtmosCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMConstAtmosCoupler::initFromOptions()\n");

  ierr = PISMAtmosphereCoupler::initFromOptions(g); CHKERRQ(ierr);
  
  if (initializeFromFile) {
    char filename[PETSC_MAX_PATH_LEN];
    LocalInterpCtx* lic;

    ierr = findPISMInputFile((char*) filename, lic); CHKERRQ(ierr); // allocates lic
    ierr = verbPrintf(2, g->com, 
       "initializing constant atmospheric climate: reading net surface mass\n"
       "  balance 'acab' and absolute surface temperature 'artm' from %s ...\n",
       filename); CHKERRQ(ierr); 

    ierr = vsurfmassflux.regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = vsurftemp.regrid(filename, *lic, true); CHKERRQ(ierr);

    delete lic;
  }

  printIfDebug("ending PISMConstAtmosCoupler::initFromOptions()\n");
  return 0;
}



/*******************  ATMOSPHERE:  PISMSnowModelAtmosCoupler ********************/

PISMSnowModelAtmosCoupler::PISMSnowModelAtmosCoupler() : PISMAtmosphereCoupler() {
  monthlysnowtemps = NULL;
  mbscheme = NULL;
  snowtempSummerWarming = config.get("-pdd_summer_warming"); // FIXME: change name to -snow_temp_summer_warming
  snowtempSummerPeakDay = config.get("-pdd_summer_peak_day"); // FIXME: change name to -snow_temp_summer_peak_day
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


bool PISMSnowModelAtmosCoupler::optionsChooseSnowModel() {
  const int N = 10;
  const char optNames[N][30] = {
           "-pdd",
           "-pdd_factor_snow",
           "-pdd_factor_ice",
           "-pdd_refreeze", 
           "-pdd_rand",
           "-pdd_rand_repeatable",
           "-pdd_std_dev",
           "-pdd_summer_warming",
           "-pdd_summer_peak_day",
           "-pdd_monthly_temps"        };
  PetscTruth check;
  for (int k = 0; k < N; k++) {
    check_option(optNames[k], check);
    if (check == PETSC_TRUE)  return true;
  }
  return false;
}


PetscErrorCode PISMSnowModelAtmosCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  PetscTruth   optSet;
  printIfDebug("entering PISMSnowModelAtmosCoupler::initFromOptions()\n");

  ierr = PISMAtmosphereCoupler::initFromOptions(g); CHKERRQ(ierr);
  
  if (!optionsChooseSnowModel()) {
    ierr = verbPrintf(1, g->com, 
       "WARNING:  PISMSnowModelAtmosCoupler::initFromOptions() called but user seems not to want it ...\n");
       CHKERRQ(ierr);
  }

  // check on whether we should read monthly snow temperatures from file  
  char monthlyTempsFile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-pdd_monthly_temps", 
             monthlyTempsFile, PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet == PETSC_TRUE) {
    ierr = verbPrintf(2,grid->com,
       "  snow temperature cycle: reading monthly temperatures from file %s ...\n",
       monthlyTempsFile); CHKERRQ(ierr);
    monthlysnowtemps = new MonthlyDataMaps;
    // puts month-by-month "reading ..." message to stdout:
    ierr = monthlysnowtemps->initFromFile(
             "monsnowtemp", // expects short names monsnowtemp1,...,monsnowtemp12
             "",            // no standard_name?; FIXME: check on this
             monthlyTempsFile, grid); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid->com,
       "  snow temperature cycle: using temperature parameterization ...\n"); CHKERRQ(ierr);
    // monthlysnowtemps == NULL from now on
  }

  // check if user wants default or random PDD; initialize mbscheme to one of these PDDs
  PetscTruth  pddRandSet, pddRepeatableSet;
  ierr = check_option("-pdd_rand", pddRandSet); CHKERRQ(ierr);
  ierr = check_option("-pdd_rand_repeatable", pddRepeatableSet); CHKERRQ(ierr);
  if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
    mbscheme = new PDDrandMassBalance((pddRepeatableSet == PETSC_TRUE));
  } else { // default case
    mbscheme = new PDDMassBalance;
  }
  ierr = mbscheme->init(); CHKERRQ(ierr);

  // check options for parameter values related to temp parameterization
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_warming", // FIXME: change name to -snow_temp_summer_warming
             &snowtempSummerWarming, &optSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_peak_day",// FIXME: change name to -snow_temp_summer_peak_day
             &snowtempSummerPeakDay, &optSet); CHKERRQ(ierr);

  // read snow precipitation rate and snow temperature
  char  filename[PETSC_MAX_PATH_LEN];
  LocalInterpCtx* lic;
  ierr = findPISMInputFile((char*)filename, lic); CHKERRQ(ierr); // allocates and initializes lic

  // snow temperatures used within mass balance model and from the temperature parameterization;
  //   these are ignored if monthly temps are used
  ierr = vsnowtemp_ma.create(*g, "snowtemp_ma", false); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.set_attrs(
            "climate_state",
            "mean annual snow temperature used in mass balance scheme",
            "K",
            ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = vsnowtemp_ma.set(273.15); CHKERRQ(ierr);  // merely a default value
  ierr = vsnowtemp_mj.create(*g, "snowtemp_mj", false); CHKERRQ(ierr);
  ierr = vsnowtemp_mj.set_attrs(
            "climate_state",
            "mean July snow temperature used in mass balance scheme",
            "K",
            ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = vsnowtemp_mj.set(273.15); CHKERRQ(ierr);  // merely a default value

  // mean annual ice equivalent snow accumulation rate (before melt, and not including rain)
  ierr = vsnowprecip.create(*g, "snowaccum", false); CHKERRQ(ierr); // FIXME: change name to snowprecip
  ierr = vsnowprecip.set_attrs(
            "climate_state", 
            "mean annual ice-equivalent snow accumulation rate",
	    "m s-1", 
	    "");  // no CF standard_name ??
	    CHKERRQ(ierr);
  ierr = vsnowprecip.set_glaciological_units("m year-1");
  vsnowprecip.write_in_glaciological_units = true;
  ierr = vsnowprecip.set(0.0); CHKERRQ(ierr);  // merely a default value

  // read snow precipitation rate from file
  ierr = verbPrintf(2, g->com, 
      "  reading mean annual ice-equivalent snow accumulation rate 'snowaccum' from %s ... \n",
      filename); CHKERRQ(ierr); 
  // this will continue with "not found" message if not found:
  ierr = vsnowprecip.regrid(filename, *lic, false); CHKERRQ(ierr);

  delete lic;

  printIfDebug("ending PISMSnowModelAtmosCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMSnowModelAtmosCoupler::writeCouplingFieldsToFile(const char *filename) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::writeCouplingFieldsToFile(filename); CHKERRQ(ierr);

  // FIXME: compute save from yearly cycle or by interpolating monthlies?
  // here we just save the mean annual and mean july only valid if a temperature parameterization was used
  ierr = vsnowtemp_ma.write(filename, NC_FLOAT); CHKERRQ(ierr);  
  ierr = vsnowtemp_mj.write(filename, NC_FLOAT); CHKERRQ(ierr);  

  ierr = vsnowprecip.write(filename, NC_FLOAT); CHKERRQ(ierr);
  return 0;
}


/*!
  See formulas (1) and (2) and Table 3 in \ref Faustoetal2009.
  
  Default values, stored in src/pism_config.cdl, use lines 
  "Best annual fit: This study with \f$\kappa_{\text{ma}}\f$" 
  and "Best July fit: This study with \f$\kappa_{\text{ma}}\f$"
  from Table 3.
 */
PetscErrorCode PISMSnowModelAtmosCoupler::parameterizationToUpdateSnowTemp(
                                IceInfoNeededByAtmosphereCoupler* info) {
  PetscErrorCode ierr;
  const PetscScalar 
    d_ma     = config.get("snow_temp_d_ma"),      // K
    gamma_ma = config.get("snow_temp_gamma_ma"),  // K m-1
    c_ma     = config.get("snow_temp_c_ma"),      // K (degN)-1
    kappa_ma = config.get("snow_temp_kappa_ma"),  // K (degW)-1
    d_mj     = config.get("snow_temp_d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = config.get("snow_temp_gamma_mj"),
    c_mj     = config.get("snow_temp_c_mj"),
    kappa_mj = config.get("snow_temp_kappa_mj");
  
  PetscScalar **lat_degN, **lon_degE, **h, **T_ma, **T_mj;
  ierr = info->surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = info->lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = info->lon->get_array(lon_degE); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);

  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      T_ma[i][j] = d_ma + gamma_ma * h[i][j] + c_ma * lat_degN[i][j] + kappa_ma * (-lon_degE[i][j]);
      T_mj[i][j] = d_mj + gamma_mj * h[i][j] + c_mj * lat_degN[i][j] + kappa_mj * (-lon_degE[i][j]);
    }
  }
  
  ierr = info->surfelev->end_access();   CHKERRQ(ierr);
  ierr = info->lat->end_access(); CHKERRQ(ierr);
  ierr = info->lon->end_access(); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  return 0;
}


//! Compute the ice surface mass flux from the snow temperature parameterization, a stored map of snow precipication rate, and a choice of mass balance scheme (typically a PDD).
/*!
First this method recomputes the snow temperature, either from a parameterization or from interpolating
stored monthly snow temperature maps.  More precisely, if there are no monthly temperature maps then
useTemperatureParameterizationToUpdateSnowTemp() is called .  In either case, at each point on the surface
a temperature time series with short (weekly or less) time steps is generated.
When there are no monthly temperatures, there is a parameterized yearly cycle 
of temperature and additional weather-related variability according to a normally distributed 
random temperature change for each week and grid point.

At each point on the ice surface a temperature time series is used by a LocalMassBalance object
to compute the number of positive degree days.  There are two such schemes, derived 
classes of LocalMassBalance.  There is a deterministic
(expectation) default method and random (monte carlo) method.  The standard deviation of
the temperature variation change can be controlled by option <tt>-pdd_std_dev</tt>.  
The deterministic method computes only the expected number of positive degree days
for that amount of variability \ref CalovGreve05; it is chosen by option <tt>-pdd</tt>.
The random method uses pseudo-random numbers to simulate the variability and then directly 
sums the number of positive degree days.  It is chosen by either <tt>-pdd_rand</tt> or 
<tt>-pdd_rand_repeatable</tt>.

A more realistic pattern for the variability of surface melting would have correlation 
with appropriate spatial and temporal ranges.

Alternatively, if option <tt>-pdd_monthly_temps</tt> provides a NetCDF file 
with the 12 monthly temperature maps, the temperature on a given day at a given
location is found by linear interpolation (in time) of the monthly temps.

The surface mass balance is computed from the stored map of snow precipitation rate
by a call to the getMassFluxFromTemperatureTimeSeries() method of the LocalMassBalance object.
 */
PetscErrorCode PISMSnowModelAtmosCoupler::updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type iceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvsmf) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
     t_years, dt_years, iceInfoNeeded, pvsmf); CHKERRQ(ierr);
 
  IceInfoNeededByAtmosphereCoupler* info = (IceInfoNeededByAtmosphereCoupler*) iceInfoNeeded;

  // Set up temperature time series.  Because Calov-Greve method uses Simpson's rule to do integral,
  // we choose the number of times to be odd.  Accuracy for that method, with a smooth yearly cycle,
  // suggests that at least 53 evals per year (i.e. approximately
  // weekly) should be sufficiently accurate.  We use the same number of times for the random PDD
  // method, too.
  PetscInt          Ntemps = (int) ceil(52 * (dt_years) + 1);
  if (Ntemps < 3) Ntemps = 3;
  if ((Ntemps % 2) == 0)  Ntemps++;  // guarantee it is odd
  // times at which we compute snow temps are t, t+dttemps, ..., t+(Ntemps-1)dttemps;
  const PetscScalar dttemps = (dt_years * secpera) / ((PetscScalar) Ntemps);
  PetscScalar       *Tseries = new PetscScalar[Ntemps];
  
  const PetscScalar radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
                    sperd = 8.64e4, // exact number of seconds per day
                    julydaysec = sperd * config.get("snow_temp_july_day");
  
  // get access to snow temperatures (and update parameterization if appropriate)
  PetscScalar **T_ma, **T_mj, **monthlytemp[12], **snowrate, **smflux;
  if (monthlysnowtemps != NULL) {
    for (PetscInt m = 0; m < 12; ++m) {
      ierr = monthlysnowtemps->vdata[m].get_array(monthlytemp[m]); CHKERRQ(ierr);
    }
  } else {
    ierr = parameterizationToUpdateSnowTemp(info); CHKERRQ(ierr);
  }
  ierr = vsnowtemp_ma.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);
  ierr = vsnowprecip.get_array(snowrate);  CHKERRQ(ierr);
  ierr = vsurfmassflux.get_array(smflux);  CHKERRQ(ierr);

  // run through grid and build temperature time series and then get surface balance at each point
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {

      if (monthlysnowtemps != NULL) {
        // FIXME: fill in method for extracting temperature time series from monthly data
      } else {
        // default case: use corrected formula (4) in \ref Faustoetal2009
        for (PetscInt k = 0; k < Ntemps; ++k) {
          const PetscScalar
            tk = (t_years - floor(t_years)) * secpera + ((PetscScalar) k) * dttemps;
          Tseries[k] = T_ma[i][j] + (T_mj[i][j] - T_ma[i][j]) * cos(radpersec * (tk - julydaysec));
        }
      }


// FIXME: NEEDS WORK FROM HERE
#if 0

      PetscScalar pdd_sum = 0.0;  // units of day^-1 K-1

      const PetscScalar summer_warming = getSummerWarming(h[i][j],lat[i][j],amstemp[i][j]);

      // use one of the two methods for computing the number of positive degree days
      // at the given i,j grid point for the duration of time step (=dt)
      if (pddRandGen != NULL) { // since random stuff set up, do monte carlo
        pdd_sum = 0.0;
        // compute # of pos deg day:
        for (PetscInt day = intstartday; day < intstartday + num_days; day++){ 
          PetscScalar mytemp;
          if (readMonthlyTempsFromFile) {
            PetscInt currMonthInd, nextMonthInd;
            PetscScalar indexday, dummy, lambda;
            indexday = 364.24 * modf(((PetscScalar) day)/365.24, &dummy); 
            ierr = getMonthIndicesFromDay(indexday, 
                       currMonthInd, nextMonthInd, lambda); CHKERRQ(ierr);
            mytemp = getTemperatureFromMonthlyData(
                       smonthtemp[currMonthInd], smonthtemp[nextMonthInd], lambda, i, j);
          } else {
            mytemp = getTemperatureFromYearlyCycle(summer_warming, amstemp[i][j], (PetscScalar) day);
          }
          const double randadd = gsl_ran_gaussian(pddRandGen, pddStdDev);
          const PetscScalar temp = mytemp + (PetscScalar) randadd;
          insttemp[i][j] = temp;
          if (temp > 273.15)   pdd_sum += (temp - 273.15);
        }
      } else { // default Calov-Greve method; apply Simpson's rule for integral
        pdd_sum = 0.0;
        for (PetscInt m = 0; m < CGsumcount; ++m) {
          // Simpson's is: (h/3) * sum([1 4 2 4 2 4 ... 4 1] .* [f(x0) f(x1) ... f(xN)])
          PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
          if ( (m == 0) || (m == (CGsumcount - 1)) )  coeff = 1.0;
          const PetscScalar day = startday + m * CGsumstepdays;
          PetscScalar temp;
          if (readMonthlyTempsFromFile) {
            PetscInt currMonthInd, nextMonthInd;
            PetscScalar indexday, dummy, lambda;
            indexday = 364.24 * modf(((PetscScalar) day)/365.24, &dummy); 
            ierr = getMonthIndicesFromDay(indexday, currMonthInd, nextMonthInd, lambda); CHKERRQ(ierr);
            temp = getTemperatureFromMonthlyData(
                       smonthtemp[currMonthInd], smonthtemp[nextMonthInd], lambda, i, j);
          } else {
            temp = getTemperatureFromYearlyCycle(summer_warming, amstemp[i][j], day);
          }
          insttemp[i][j] = temp;
          pdd_sum += coeff * CalovGreveIntegrand(temp);  // temp in K
        }
        pdd_sum = (CGsumstepdays / 3.0) * pdd_sum;
      }

      // now that we have the number of PDDs, compute mass balance from snow rate
      smflux[i][j] = getSurfaceBalanceFromSnowAndPDD(
                            saccum[i][j], dt_years * secpera, pdd_sum);

#endif

    }
  }
  
  if (monthlysnowtemps != NULL) {
    for (PetscInt m = 0; m < 12; ++m) {
      ierr = monthlysnowtemps->vdata[m].end_access(); CHKERRQ(ierr);
    }
  }
  ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);
  ierr = vsnowprecip.end_access();  CHKERRQ(ierr);
  ierr = vsurfmassflux.end_access(); CHKERRQ(ierr);

  delete [] Tseries;
  return 0;
}


/*******************  ATMOSPHERE:  PISMMonthlyTempsAtmosCoupler ********************/

PISMMonthlyTempsAtmosCoupler::PISMMonthlyTempsAtmosCoupler() : PISMAtmosphereCoupler() {
  readMonthlyTempsFromFile = false;
  strcpy(monthlyTempsFile,""); // zero length file name; causes error if readMonthlyTemps() called
}


PISMMonthlyTempsAtmosCoupler::~PISMMonthlyTempsAtmosCoupler() {
  for (PetscInt j = 0; j < 12; ++j) {
    vmonthlytemp[j].destroy();  // destroys if allocated (created)
  }
}


//! Initializes by reading monthly temperatures from the PISM input file.
/*!
Stored temperatures must have names \c monsnowtemp1, ...,\c monsnowtemp12 and be in units of K.
Call setMonthlyTempsFilename() and make sure readMonthlyTempsFromFile == true before 
using this initFromOptions().
 */
PetscErrorCode PISMMonthlyTempsAtmosCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMMonthlyTempsAtmosCoupler::initFromOptions()\n");

  ierr = PISMAtmosphereCoupler::initFromOptions(g); CHKERRQ(ierr); // sets grid and metadata

  if (readMonthlyTempsFromFile) {
    ierr = readMonthlyTemps(); CHKERRQ(ierr);
  }

  printIfDebug("ending PISMMonthlyTempsAtmosCoupler::initFromOptions()\n");
  return 0;
}


//! Write monthly temperatures to a prepared file.
/*!
Adds \c monsnowtemp1, ...,\c monsnowtemp12 after other PISMAtmosphereCoupler fields.
 */
PetscErrorCode PISMMonthlyTempsAtmosCoupler::writeCouplingFieldsToFile(const char *filename) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::writeCouplingFieldsToFile(filename); CHKERRQ(ierr);

  for (PetscInt j = 0; j < 12; ++j) {
    if (vmonthlytemp[j].was_created()) {
      ierr = vmonthlytemp[j].write(filename, NC_FLOAT); CHKERRQ(ierr);
    }
  }
  return 0;
}


//! Set the name of the NetCDF file from which we read monthly temperatures.
PetscErrorCode PISMMonthlyTempsAtmosCoupler::setMonthlyTempsFilename(const char* filename) {
  strcpy(monthlyTempsFile,filename);
  return 0;
}


//! Read monthly temperatures from a prepared file.
/*!
Reads \c monsnowtemp1, ...,\c monsnowtemp12 from file with name monthlyTempsFile.
 */
PetscErrorCode PISMMonthlyTempsAtmosCoupler::readMonthlyTemps() {
  PetscErrorCode ierr;

  if (!readMonthlyTempsFromFile) {
    ierr = PetscPrintf(grid->com, 
       "PISMMonthlyTempsAtmosCoupler ERROR:  readMonthlyTempsFromFile == false\n"); CHKERRQ(ierr);
    PetscEnd();
  }
  if (strlen(monthlyTempsFile) == 0) {
    ierr = PetscPrintf(grid->com, 
       "PISMMonthlyTempsAtmosCoupler ERROR:  empty filename for file from which to read\n"
       "                                     temps (monthlyTempsFile)\n"); CHKERRQ(ierr);
    PetscEnd();
  }
  
  // find file and set up info so regrid works
  NCTool nc(grid);
  grid_info gi;
  ierr = nc.open_for_reading(monthlyTempsFile); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  LocalInterpCtx lic(gi, NULL, NULL, *grid); // 2D only; destructed by end of scope

  // for each month, create an IceModelVec2 and assign attributes
  for (PetscInt j = 0; j < 12; ++j) {
    char monthlyTempName[20], mTstring[100];
    snprintf(monthlyTempName, 20, "monsnowtemp%d", j+1);
    ierr = vmonthlytemp[j].create(*grid, monthlyTempName, false); // global; no ghosts
       CHKERRQ(ierr);
    snprintf(mTstring, 100, 
             "temperature at ice surface during month %d of {1,..,12}", j+1);
    ierr = vmonthlytemp[j].set_attrs(
               "climate_state",
               mTstring, // note simplified, not-very-specific long name
               "K",
               ""); // CF standard name?  may exist when derived class has additional semantics
               CHKERRQ(ierr);
    ierr = verbPrintf(2, grid->com, 
       "  reading month %d surface temperature '%s' ...\n", j+1, monthlyTempName); CHKERRQ(ierr); 
    ierr = vmonthlytemp[j].regrid(monthlyTempsFile, lic, true); CHKERRQ(ierr); // it *is* critical
  }
  return 0;
}


/*!
Used for indexing the monthly surface temperature data.
Returns indices in {0,..,11} for the current and next months.  Returns linear interpolation
coefficient lambda (with 0 <= lambda < 1) used by getTemperatureFromMonthlyData().
 */
PetscErrorCode PISMMonthlyTempsAtmosCoupler::getMonthIndicesFromDay(
       PetscScalar day, PetscInt &curr, PetscInt &next, PetscScalar &lambda) {
  PetscErrorCode ierr;
  if ((day < 0) || (day > 365.24)) {
    ierr = PetscPrintf(grid->com, 
       "PISMMonthlyTempsAtmosCoupler ERROR:  invalid day = %.5f in getMonthIndicesFromDay();\n"
       "      should be in range 0 <= day <= 365.24\n",day); CHKERRQ(ierr);
    PetscEnd();
  }
  PetscScalar month = 12.0 * day / 365.24;  // 0 <= month < 12
  lambda = month - floor(month);            // 0 <= lambda < 1
  curr = (int) floor(month);                // curr in {0,1,2,...,10,11}
  next = curr + 1;
  if (next == 12)   next = 0;               // next in {0,1,2,...,10,11}
  return 0;
}


/*!
Linearly interpolates between stored monthly temps.  Call getMonthIndicesFromDay()
first to get temporal index into monthly data and for linear interpolation factor lambda.
 */
PetscScalar PISMMonthlyTempsAtmosCoupler::getTemperatureFromMonthlyData(
       PetscScalar **currMonthTemps, PetscScalar **nextMonthTemps, PetscScalar lambda,
       PetscInt i, PetscInt j) {
  return (1.0 - lambda) * currMonthTemps[i][j] + lambda * nextMonthTemps[i][j]
       + TsOffset;
}


/*******************  OCEAN:  PISMOceanCoupler ********************/

PISMOceanCoupler::PISMOceanCoupler() : PISMClimateCoupler() {
  reportInitializationToStdOut = true;  // derived classes etc. can turn off before calling
                                        // initFromOptions(), but its on by default
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
PetscErrorCode PISMOceanCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMOceanCoupler::initFromOptions()\n");

  ierr = PISMClimateCoupler::initFromOptions(g); CHKERRQ(ierr);

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
    dSLforcing = new PISMClimateForcing;
    int stat, ncid = 0;
    ierr = verbPrintf(2, grid->com, 
         "  reading delta sea level data from forcing file %s ...\n", 
         dSLfile); CHKERRQ(ierr);
    if (grid->rank == 0) {    stat = nc_open(dSLfile, 0, &ncid); CHKERRQ(nc_check(stat));   }
    ierr = dSLforcing->readClimateForcingData(grid->com, grid->rank, 
             ncid, grid->year,PCF_DELTA_SEA_LEVEL); CHKERRQ(ierr);
    if (grid->rank == 0) {    stat = nc_close(ncid); CHKERRQ(nc_check(stat));   }
  }

  printIfDebug("ending PISMOceanCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMOceanCoupler::writeCouplingFieldsToFile(const char *filename) {
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
  const PetscScalar /*t_years*/, const PetscScalar /*dt_years*/, 
  void */*iceInfoNeeded*/, IceModelVec2* &pvsbmf) {

  if (vshelfbasemassflux.was_created())
    pvsbmf = &vshelfbasemassflux;
  else {  
    SETERRQ(1,"vvshelfbasemassflux not created in PISMOceanCoupler::updatehelfBaseMassFluxAndProvide()");
  }

  return 0;
}


//! Provides access to vshelfbasetemp.  No update of vshelfbasetemp.  Real ocean models in derived classes will update.
PetscErrorCode PISMOceanCoupler::updateShelfBaseTempAndProvide(
                  const PetscScalar /*t_years*/, const PetscScalar /*dt_years*/, 
                  void */*iceInfoNeeded*/, IceModelVec2* &pvsbt) {
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
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded) {
  PetscErrorCode ierr;
  IceModelVec2* ignored;
  ierr = updateShelfBaseMassFluxAndProvide(t_years, dt_years, iceInfoNeeded, ignored); CHKERRQ(ierr);
  ierr = updateShelfBaseTempAndProvide(t_years, dt_years, iceInfoNeeded, ignored); CHKERRQ(ierr);
  ierr = updateSeaLevelElevation(t_years, dt_years, NULL); CHKERRQ(ierr);
  return 0;
}


//! Updates the sea level (using -dSLforcing, if it is on) and sets \c new_sea_level (if not NULL).
PetscErrorCode PISMOceanCoupler::updateSeaLevelElevation(PetscReal t_years, PetscReal /*dt_years*/,
							 PetscReal *new_sea_level) {
  PetscErrorCode ierr;

  if (dSLforcing != PETSC_NULL) {
    // read the new sea level (delta from modern)
    ierr = dSLforcing->updateFromClimateForcingData(t_years, &seaLevel); CHKERRQ(ierr);
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


PetscErrorCode PISMConstOceanCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMConstOceanCoupler::initFromOptions()\n");

  ierr = PISMOceanCoupler::initFromOptions(g); CHKERRQ(ierr);
  
  if (reportInitializationToStdOut) {
    ierr = verbPrintf(2, g->com, 
       "initializing constant sub-ice shelf ocean climate: heat flux from ocean\n"
       "  set to %.3f W m-2; this determines mass balance; ice shelf base temperature\n"
       "  set to pressure-melting temperature ...\n",
       constOceanHeatFlux); CHKERRQ(ierr); 
  }

  printIfDebug("ending PISMConstOceanCoupler::initFromOptions()\n");
  return 0;
}


//! By default, does not write vshelfbasetemp and vshelfbasemassflux fields.
PetscErrorCode PISMConstOceanCoupler::writeCouplingFieldsToFile(const char *filename) {
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
                  const PetscScalar /*t_years*/, const PetscScalar /*dt_years*/, 
                  void *iceInfoNeeded, IceModelVec2* &pvsbt) {
  PetscErrorCode ierr;

  // ignores everything from IceModel except ice thickness; also ignores t_years and dt_years
  const PetscScalar icerho       = 910.0,   // kg/m^3   ice shelf mean density
                    oceanrho     = 1028.0,  // kg/m^3   sea water mean density
                    beta_CC_grad = 8.66e-4, // K/m      Clausius-Clapeyron gradient
                    T0           = 273.15;  // K        triple point = naively assumed
                                            //          sea water temperature at sea level

  IceInfoNeededByOceanCoupler* info = (IceInfoNeededByOceanCoupler*) iceInfoNeeded;
  PetscScalar **H, **temp;

  ierr = info->thk->get_array(H);   CHKERRQ(ierr);
  ierr = vshelfbasetemp.get_array (temp); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      const PetscScalar shelfbaseelev = - (icerho / oceanrho) * H[i][j];
      // temp is set to melting point at depth
      temp[i][j] = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so this is below T0
    }
  }
  ierr = info->thk->end_access(); CHKERRQ(ierr);
  ierr = vshelfbasetemp.end_access(); CHKERRQ(ierr);
  
  pvsbt = &vshelfbasetemp;
  return 0;                                 
}


PetscErrorCode PISMConstOceanCoupler::updateShelfBaseMassFluxAndProvide(
                  const PetscScalar /*t_years*/, const PetscScalar /*dt_years*/, 
                  void */*iceInfoNeeded*/, IceModelVec2* &pvsbmf) {
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

