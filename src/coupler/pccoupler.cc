// Copyright (C) 2009 Ed Bueler, Ricarda Winkelmann and Constantine Khroulev
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
#include "../base/iceModel.hh"
#include "../base/iceModelVec.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/nc_util.hh"
#include "pccoupler.hh"


/******************* VIRTUAL BASE CLASS:  PISMClimateCoupler ********************/

PISMClimateCoupler::PISMClimateCoupler() {
  grid = NULL;
}


PISMClimateCoupler::~PISMClimateCoupler() { 
}


PetscErrorCode PISMClimateCoupler::initFromOptions(IceGrid* g) {
  grid = g;
  return 0;
}


PetscErrorCode PISMClimateCoupler::findPISMInputFile(char* filename, LocalInterpCtx* &lic) {
  if (grid == NULL) {
    SETERRQ(1,"findPISMInputFile(): grid not initialized");
  }

  PetscErrorCode ierr;
  PetscTruth ifSet, bifSet;
  PetscTruth i_set, boot_from_set;

  // OLD OPTIONS
  ierr = PetscOptionsHasName(PETSC_NULL, "-if", &ifSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-bif", &bifSet); CHKERRQ(ierr);
  // NEW OPTIONS
  ierr = PetscOptionsHasName(PETSC_NULL, "-i", &i_set); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-boot_from", &boot_from_set); CHKERRQ(ierr);
  // warnings to let users get used to the change:
  if (ifSet) {
    ierr = verbPrintf(2, grid->com,
		      "PISM WARNING: '-if' command line option is deprecated. Please use '-i' instead.\n");
    CHKERRQ(ierr);
  }
  if (bifSet) {
    ierr = verbPrintf(2, grid->com, 
		      "PISM WARNING: '-bif' command line option is deprecated. Please use '-boot_from' instead.\n");
    CHKERRQ(ierr);
  }

  char i_file[PETSC_MAX_PATH_LEN], boot_from_file[PETSC_MAX_PATH_LEN];
  if (i_set) {
    if (boot_from_set) {
      ierr = PetscPrintf(grid->com,
			 "PISM ERROR: both '-i' and '-boot_from' are used. Exiting...\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    if (ifSet) {		// OLD OPTION
      ierr = PetscPrintf(grid->com,
			 "PISM ERROR: both '-i' and '-if' are used. Ignoring '-if'...\n"); CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetString(PETSC_NULL, "-i", i_file, 
				 PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
    strcpy(filename, i_file);
  }
  else if (boot_from_set) {
    if (ifSet) {		// OLD OPTION
      ierr = PetscPrintf(grid->com,
			 "PISM ERROR: both '-if' and '-boot_from' are used. Exiting...\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    
    ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from", boot_from_file, 
				 PETSC_MAX_PATH_LEN, &boot_from_set); CHKERRQ(ierr);
    strcpy(filename, boot_from_file);
  }
  else if (ifSet) {		// OLD OPTION
    if (bifSet) {
      ierr = PetscPrintf(grid->com,
			 "PISM ERROR: both '-bif' and '-if' are used. Exiting...\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    ierr = PetscOptionsGetString(PETSC_NULL, "-if", i_file, 
				 PETSC_MAX_PATH_LEN, &ifSet); CHKERRQ(ierr);
    strcpy(filename, i_file);
  }
  else if (bifSet) {		// OLD OPTION
    ierr = PetscOptionsGetString(PETSC_NULL, "-bif", boot_from_file, 
				 PETSC_MAX_PATH_LEN, &bifSet); CHKERRQ(ierr);
    strcpy(filename, boot_from_file);
  } else {
    PetscPrintf(grid->com, "PISM ERROR: no -i and no -boot_from specified. Exiting...\n");
    PetscEnd();
  }

  // file name now contains name of PISM input file;
  //   now we check it is really there and, if so, we read the dimensions of
  //   the computational grid so that we can set up a LocalInterpCtx for further
  //   reading of actual climate data from the file
  bool file_exists = false;
  NCTool nc(grid);
  grid_info gi;
  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid->com,
             "PISMClimateCoupler ERROR: Can't open file '%s'.\n",
             filename); CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // *caller* of findPISMInputFile() is in charge of destroying
  lic = new LocalInterpCtx(gi, NULL, NULL, *grid); // 2D only

  return 0;
}


PetscErrorCode PISMClimateCoupler::updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, void *iceInfoNeeded) {
  SETERRQ(1,"PCC ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
  return 0;
}


PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(const char *filename) {
  SETERRQ(1,"PCC ERROR:  this method is VIRTUAL in PISMClimateCoupler and not implemented");
  return 0;
}


/******************* ATMOSPHERE:  PISMAtmosphereCoupler ********************/

PISMAtmosphereCoupler::PISMAtmosphereCoupler() : PISMClimateCoupler() {
}


PISMAtmosphereCoupler::~PISMAtmosphereCoupler() {
  vsurfmassflux.destroy();
  vsurftemp.destroy();
}


/*!
Derived class implementations will check user options to configure the PISMAtmosphereCoupler.
This version allocates space and sets attributes for the two essential fields.
 */
PetscErrorCode PISMAtmosphereCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;

  ierr = PISMClimateCoupler::initFromOptions(g); CHKERRQ(ierr);
  
  // short names "acab" and "artm" match GLIMMER (& CISM, presumably)
  
  // mean annual net ice equivalent surface mass balance rate
  //   FIXME re create(): IceModelVec2 is local even though ghosts are not needed;
  //   this is so that vsurfmassflux duplicates IceModel::vAccum; later we
  //   can make global with no ghosts
  ierr = vsurfmassflux.create(*g, "acab", true); CHKERRQ(ierr);
  ierr = vsurfmassflux.set_attrs(
            "climate_state", 
            "mean annual net ice equivalent accumulation (ablation) rate",
	    "m s-1", 
	    "land_ice_surface_specific_mass_balance");  // CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfmassflux.set_glaciological_units("m year-1");
  vsurfmassflux.write_in_glaciological_units = true;
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr);  // merely a default value

  // annual mean air temperature at "ice surface", at level below all firn processes
  // possibly should be reported in deg C; would require shift version of glaciological_units
  //   FIXME re create(): ditto above; vsurftemp follows vTs
  ierr = vsurftemp.create(*g, "artm", true); CHKERRQ(ierr);
  ierr = vsurftemp.set_attrs(
            "climate_state",
            "temperature at ice surface but below firn",
            "K", 
            NULL);  // PROPOSED CF standard_name = land_ice_temperature_below_firn
            CHKERRQ(ierr);
  ierr = vsurftemp.set(273.15); CHKERRQ(ierr);  // merely a default value
  return 0;
}


PetscErrorCode PISMAtmosphereCoupler::writeCouplingFieldsToFile(const char *filename) {
  PetscErrorCode ierr;
  
  // We assume file is prepared in the sense that it exists and that global attributes 
  //   are already written.  See IceModel::dumpToFile() for how main PISM output file is
  //   prepared.  Note calls below do handle opening and closing the file.
  //   Note we write in FLOAT not DOUBLE because these are expected to be for diagnosis,
  //   not (e.g.) restart etc.
  ierr = vsurfmassflux.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = vsurftemp.write(filename, NC_FLOAT); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsmf) {
  SETERRQ(1,"VIRTUAL in PISMAtmosphereCoupler ... not implemented");
  return 0;
}


PetscErrorCode PISMAtmosphereCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvst) {
  SETERRQ(1,"VIRTUAL in PISMAtmosphereCoupler ... not implemented");
  return 0;
}


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


PetscErrorCode PISMConstAtmosCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;

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
  } else {
    ierr = verbPrintf(2, g->com, 
       "constant atmospheric climate already initialized; not read from any file ...\n");
       CHKERRQ(ierr); 
  }
  return 0;
}


//! Just provides access.  Nothing to update.
PetscErrorCode PISMConstAtmosCoupler::updateSurfMassFluxAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvsmf) {
  pvsmf = &vsurfmassflux;
  return 0;
}


//! Just provides access.  Nothing to update
PetscErrorCode PISMConstAtmosCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvst) {
  pvst = &vsurftemp;
  return 0;
}


/*******************  OCEAN:  PISMOceanCoupler ********************/

PISMOceanCoupler::PISMOceanCoupler() : PISMClimateCoupler() {
}

