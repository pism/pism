// Copyright (C) 2009 Ed Bueler and Ricarda Winkelmann
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


PetscErrorCode PISMClimateCoupler::findPISMInputFile(
                     char* *filename, LocalInterpCtx* &lic) {
  if (grid == NULL) {
    SETERRQ(1,"findPISMInputFile(): grid not initialized");
  }

  PetscErrorCode ierr;
  PetscTruth ifSet, bifSet;
  char if_file[PETSC_MAX_PATH_LEN], bif_file[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", if_file, 
                               PETSC_MAX_PATH_LEN, &ifSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", bif_file, 
                               PETSC_MAX_PATH_LEN, &bifSet); CHKERRQ(ierr);
  if (bifSet == PETSC_TRUE) {
    strcpy(*filename,bif_file);
  } else if (ifSet == PETSC_TRUE) {
    strcpy(*filename,if_file);
  } else {
    SETERRQ(2,"findPISMInputFile(): no -if and no -bif?; how did I get here?");
  }

  bool file_exists = false;
  NCTool nc(grid);
  grid_info gi;

  ierr = nc.open_for_reading(*filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid->com,
             "PISMClimateCoupler ERROR: Can't open file '%s'.\n",
             filename); CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  lic = new LocalInterpCtx(gi, NULL, NULL, *grid); // 2D only; 
  // caller is in charge of destroying
  return 0;
}


PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(const char *filename) {
  SETERRQ(1,"not implemented");
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
  ierr = vsurfmassflux.set_glaciological_units("m year-1", secpera);
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


//! Just provides access.  Generally, the surface mass flux is updated here, by atmosphere model.
PetscErrorCode PISMAtmosphereCoupler::updateSurfMassFluxAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  IceModelVec2 mask, IceModelVec2 surface_elev,
                  IceModelVec2* &pvsmf) {
  pvsmf = &vsurfmassflux;
  return 0;
}


//! Just provides access.  Generally, the surface temp is updated here, by atmosphere model.
PetscErrorCode PISMAtmosphereCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  IceModelVec2 vmask, IceModelVec2 vsurfelev,
                  IceModelVec2* &pvst) {
  pvst = &vsurftemp;
  return 0;
}


/*******************  ATMOSPHERE:  PISMConstAtmosCoupler ********************/

PISMConstAtmosCoupler::PISMConstAtmosCoupler() : PISMAtmosphereCoupler() {
}


PetscErrorCode PISMConstAtmosCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  char actualfilename[PETSC_MAX_PATH_LEN];
  char* filename=&(actualfilename[0]);
  LocalInterpCtx* lic;

  ierr = PISMAtmosphereCoupler::initFromOptions(g); CHKERRQ(ierr);

  ierr = findPISMInputFile(&filename, lic); CHKERRQ(ierr); // allocates lic

  ierr = verbPrintf(2, g->com, 
     "initializing constant atmospheric climate: reading net surface mass balance\n"
     "  'acab' and absolute surface temperature 'artm' from %s ...\n",
     filename); CHKERRQ(ierr); 

  ierr = vsurfmassflux.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical
  ierr = vsurftemp.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical

  delete lic;
  return 0;
}


/*******************  OCEAN:  PISMOceanCoupler ********************/

PISMOceanCoupler::PISMOceanCoupler() : PISMClimateCoupler() {
}

