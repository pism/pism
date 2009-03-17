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
#include "pccoupler.hh"
// note we do NOT depend on IceModel.hh; this is deliberate!


/******************* VIRTUAL BASE CLASS:  PISMClimateCoupler ********************/

PISMClimateCoupler::PISMClimateCoupler() {
  grid = NULL;
  PCCDEBUG = false;  // set to true and recompile if entry and exit messages for initFromOptions(),
                     // both base class and derived classes, are needed for debugging
}


PISMClimateCoupler::~PISMClimateCoupler() { 
}


/*!
Just set grid to provided IceGrid*.
 */
PetscErrorCode PISMClimateCoupler::initFromOptions(IceGrid* g) {
  printIfDebug("entering PISMClimateCoupler::initFromOptions()\n");
  grid = g;
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
  PetscTruth ifSet, bifSet;
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
             const PetscScalar t_years, const PetscScalar dt_years, void *iceInfoNeeded) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
  return 0;
}


//! A virtual method which writes fields associated to the derived class.
PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(const char *filename) {
  SETERRQ(1,"PISMClimateCoupler ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
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
    TsOffset = 0.0;
    delete dTforcing; // calls destructor for this PISMClimateForcing instance
    dTforcing = PETSC_NULL;
  }
}


//! Initialize a PISMAtmosphereCoupler by allocating space for surface mass flux and surface temperature variables.
/*!
Allocates space and sets attributes, including CF standard_name, for the two essential fields.
Derived class implementations will check user options to configure.

g->year must be valid before this can be called.
 */
PetscErrorCode PISMAtmosphereCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug("entering PISMAtmosphereCoupler::initFromOptions()\n");

  ierr = PISMClimateCoupler::initFromOptions(g); CHKERRQ(ierr);
  
  // short names "acab" and "artm" match GLIMMER (& CISM, presumably)
  // mean annual net ice equivalent surface mass balance rate
  ierr = vsurfmassflux.create(*g, "acab", false); CHKERRQ(ierr);
  ierr = vsurfmassflux.set_attrs(
            "climate_state", 
            "instantaneous net ice equivalent accumulation (ablation) rate",
	    "m s-1", 
	    "land_ice_surface_specific_mass_balance");  // CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfmassflux.set_glaciological_units("m year-1");
  vsurfmassflux.write_in_glaciological_units = true;
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr);  // merely a default value

  // annual mean air temperature at "ice surface", at level below all firn processes
  // possibly should be reported in deg C; would require shift version of glaciological_units
  ierr = vsurftemp.create(*g, "artm", false); CHKERRQ(ierr);
  ierr = vsurftemp.set_attrs(
            "climate_state",
            "temperature at ice surface but below firn processes",
            "K", 
            NULL);  // PROPOSED CF standard_name = land_ice_surface_temperature_below_firn
            CHKERRQ(ierr);
  ierr = vsurftemp.set(273.15); CHKERRQ(ierr);  // merely a default value

  // check user option -dTforcing for a surface temperature forcing data set
  char dTfile[PETSC_MAX_PATH_LEN];
  PetscTruth dTforceSet;
  if (dTforcing != PETSC_NULL) {
    SETERRQ(1, "dTforcing should be PETSC_NULL at start of PISMAtmosphereCoupler::initFromOptions()\n");
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
  ierr = vsurftemp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  // FIXME:  it would be good to also write the -dTforcing value to a t-dependent variable:
  //   surftempoffset(t)
  return 0;
}


//! Provides access to vsurfmassflux.  No update of vsurfmassflux.  Real atmosphere models in derived classes will update.
PetscErrorCode PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsmf) {
  if (vsurfmassflux.was_created())
    pvsmf = &vsurfmassflux;
  else {  SETERRQ(1,"vsurfmassflux not created in PISMAtmosphereCoupler::updateSurfMassFluxAndProvide()");  }
  return 0;
}


//! Updates forcing and provides access to vsurftemp.  No update of vsurftemp.  Real atmosphere models, derived classes, will update.
PetscErrorCode PISMAtmosphereCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvst) {
  PetscErrorCode ierr;

  if (dTforcing != PETSC_NULL) {
    ierr = vsurftemp.shift(-TsOffset); CHKERRQ(ierr); // return to unshifted state
    ierr = dTforcing->updateFromClimateForcingData(grid->year,&TsOffset); CHKERRQ(ierr); // read a new offset
    ierr = verbPrintf(5,grid->com,
       "PISMAtmosphereCoupler says: read TsOffset=%.6f from -dTforcing data\n",
       TsOffset); CHKERRQ(ierr);
    ierr = vsurftemp.shift(TsOffset); CHKERRQ(ierr);  // apply the offset
  }

  if (vsurftemp.was_created())
    pvst = &vsurftemp;
  else {  SETERRQ(1,"vsurftemp not created in PISMAtmosphereCoupler::updateSurfTempAndProvide()");  }

  return 0;
}


//! Calls updateSurfMassFluxAndProvide() and updateSurfTempAndProvide(), but ignors returned pointers.
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
               NULL); // CF standard name?  may exist when derived class has additional semantics
               CHKERRQ(ierr);
    ierr = verbPrintf(2, grid->com, 
       "  reading month %d surface temperature '%s' ...\n", j, monthlyTempName); CHKERRQ(ierr); 
    ierr = vmonthlytemp[j].regrid(monthlyTempsFile, lic, true); CHKERRQ(ierr); // it *is* critical
  }
  return 0;
}


/*!
Returns indices in {0,..,11} for the current and next months.  Used for indexing
the monthly surface temperature data.
 */
PetscErrorCode PISMMonthlyTempsAtmosCoupler::getMonthIndicesFromDay(
       const PetscScalar day, PetscInt &curr, PetscInt &next) {
  PetscScalar month = 12.0 * day / 365.24;
  month = month - static_cast<PetscScalar> (((int) floor(month)) % 12);
  curr = (int) floor(month);
  curr = curr % 12;
  next = curr + 1;
  if (next == 12)   next = 0;
  return 0;
}


/*!
Linearly interpolates between stored monthly temps.
 */
PetscScalar PISMMonthlyTempsAtmosCoupler::getTemperatureFromMonthlyData(
       PetscScalar **currMonthTemps, PetscScalar **nextMonthTemps,
       const PetscInt i, const PetscInt j, const PetscScalar day) {
  PetscScalar month = 12.0 * day / 365.24;
  month = month - static_cast<PetscScalar> (((int) floor(month)) % 12);
  PetscInt  curr = (int) floor(month);
  curr = curr % 12;
  return currMonthTemps[i][j] 
       + (month - (PetscScalar)curr) * nextMonthTemps[i][j];
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
	   "K", NULL); CHKERRQ(ierr);
	   // check whether a standard_name is on CF table
	   //   http://cf-pcmdi.llnl.gov/documents/cf-standard-names/current/cf-standard-name-table.html
  ierr = vshelfbasetemp.set(273.15); CHKERRQ(ierr);  // merely a default value to clear nonsense

  // ice mass balance rate at the base of the ice shelf
  ierr = vshelfbasemassflux.create(*g, "shelfbmassflux", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = vshelfbasemassflux.set_attrs(
           "climate_state", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
	   "m s-1", NULL); CHKERRQ(ierr); 
	   // check whether standard_name is on CF table
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
  // FIXME:  it would be good to also write the forcing value to a t-dependent variable:
  //   sealevel(t)
  return 0;
}


//! Updates seaLevel from forcing and provides access to vshelfbasemassflux.  No update of vshelfbasemassflux.  Real ocean models in derived classes will update.
PetscErrorCode PISMOceanCoupler::updateShelfBaseMassFluxAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsbmf) {
  PetscErrorCode ierr;

  if (dSLforcing != PETSC_NULL) {
    // read new sea level (delta from modern)
    ierr = dSLforcing->updateFromClimateForcingData(grid->year,&seaLevel); CHKERRQ(ierr);
    ierr = verbPrintf(5,grid->com,"read newSeaLevel=%.6f from -dSLforcing climate data\n",
       seaLevel); CHKERRQ(ierr);
    // comment: IceModel::updateSurfaceElevationAndMask() needs to be called before effect of changed
    //   sea level is seen in ice dynamics (e.g. on grounding line)
  }

  if (vshelfbasemassflux.was_created())
    pvsbmf = &vshelfbasemassflux;
  else {  
    SETERRQ(1,"vvshelfbasemassflux not created in PISMOceanCoupler::updatehelfBaseMassFluxAndProvide()");
  }

  return 0;
}


//! Updates seaLevel from forcing and provides access to vshelfbasetemp.  No update of vshelfbasetemp.  Real ocean models in derived classes will update.
PetscErrorCode PISMOceanCoupler::updateShelfBaseTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsbt) {
  PetscErrorCode ierr;
  //  printIfDebug("entering PISMOceanCoupler::updateShelfBaseTempAndProvide()\n");
  
  if (dSLforcing != PETSC_NULL) {
    // read new sea level (delta from modern)
    ierr = dSLforcing->updateFromClimateForcingData(grid->year,&seaLevel); CHKERRQ(ierr);
    ierr = verbPrintf(5,grid->com,"read newSeaLevel=%.6f from -dSLforcing climate data\n",
       seaLevel); CHKERRQ(ierr);
    // comment: IceModel::updateSurfaceElevationAndMask() needs to be called before effect of changed
    //   sea level is seen in ice dynamics (e.g. on grounding line)
  }

  if (vshelfbasetemp.was_created())
    pvsbt = &vshelfbasetemp;
  else {  
    SETERRQ(1,"vvshelfbasetemp not created in PISMOceanCoupler::updatehelfBaseTempAndProvide()");
  }

  //  printIfDebug("leaving PISMOceanCoupler::updateShelfBaseTempAndProvide()\n");
  return 0;
}


PetscErrorCode PISMOceanCoupler::updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded) {
  PetscErrorCode ierr;
  IceModelVec2* ignored;
  ierr = updateShelfBaseMassFluxAndProvide(t_years, dt_years, iceInfoNeeded, ignored); CHKERRQ(ierr);
  ierr = updateShelfBaseTempAndProvide(t_years, dt_years, iceInfoNeeded, ignored); CHKERRQ(ierr);
  return 0;
}


PetscReal PISMOceanCoupler::reportSeaLevelElevation() {
  return seaLevel;
}


/*******************  OCEAN:  PISMConstOceanCoupler ********************/

PISMConstOceanCoupler::PISMConstOceanCoupler() : PISMOceanCoupler() {
  constOceanHeatFlux = 0.5;   // W m-2 = J m-2 s-1; naively chosen default value
        // presumably irrelevant:  about 4 times more heating than peak of 
        //   Shapiro & Ritzwoller (2004) geothermal fluxes for Antarctica of about 130 mW/m^2
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
       "  set to %.3f W m-2 determines mass balance; ice shelf base temperature set to\n"
       "  pressure-melting temperature ...\n",
       constOceanHeatFlux); CHKERRQ(ierr); 
  }

  printIfDebug("ending PISMConstOceanCoupler::initFromOptions()\n");
  return 0;
}


PetscErrorCode PISMConstOceanCoupler::updateShelfBaseTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsbt) {
  PetscErrorCode ierr;

  // call base class to update seaLevel (if -dSLforcing); pvsbt ignored
  ierr = PISMOceanCoupler::updateShelfBaseTempAndProvide(t_years, dt_years, iceInfoNeeded, pvsbt);
            CHKERRQ(ierr);

  // ignors everything from IceModel except ice thickness; also ignors t_years and dt_years
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
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsbmf) {
  PetscErrorCode ierr;

  // call base class to update seaLevel (if -dSLforcing); pvsbmf ignored
  ierr = PISMOceanCoupler::updateShelfBaseMassFluxAndProvide(t_years, dt_years, iceInfoNeeded, pvsbmf); CHKERRQ(ierr);

  const PetscScalar icelatentheat = 3.35e5,   // J kg-1   ice latent heat capacity
                    icerho        = 910.0,    // kg m-3   ice shelf mean density
                    // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
                    meltrate      = constOceanHeatFlux / (icelatentheat * icerho); // m s-1

  // vshelfbasemassflux is positive if ice is freezing on; here it is always negative
  ierr = vshelfbasemassflux.set(- meltrate); CHKERRQ(ierr);

  pvsbmf = &vshelfbasemassflux;
  return 0;
}

