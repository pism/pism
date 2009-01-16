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
#include "pPDDcoupler.hh"


PISMPDDCoupler::PISMPDDCoupler() : PISMAtmosphereCoupler() {
}


PISMPDDCoupler::~PISMPDDCoupler() {
  vsurfaccum.destroy();
}


PetscErrorCode PISMPDDCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  char actualfilename[PETSC_MAX_PATH_LEN];
  char* filename=&(actualfilename[0]);
  LocalInterpCtx* lic;

  ierr = PISMAtmosphereCoupler::initFromOptions(g); CHKERRQ(ierr);

  // mean annual ice equivalent accumulation rate
  ierr = vsurfaccum.create(*g, "acab", true); CHKERRQ(ierr);
  ierr = vsurfaccum.set_attrs(
            "climate_state", 
            "mean annual ice equivalent accumulation rate",
	    "m s-1", 
	    NULL);  // no CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfaccum.set_glaciological_units("m year-1", secpera);
  ierr = vsurfaccum.set(0.0); CHKERRQ(ierr);  // merely a default value

  // now read two fields
  ierr = findPISMInputFile(&filename, lic); CHKERRQ(ierr);

  ierr = verbPrintf(2, g->com, 
     "initializing constant atmospheric climate and PDD: reading surface\n"
     "  accumulation 'acab' and surface temperature 'artm' from %s ... \n",
     filename); CHKERRQ(ierr); 

  ierr = vsurfaccum.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical
  ierr = vsurftemp.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical

  // now that they are read, reset vsurfaccum name to output name
  ierr = vsurfaccum.set_name("accum"); CHKERRQ(ierr);
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr); // initialize to have some values
  return 0;
}


PetscErrorCode PISMPDDCoupler::writeCouplingFieldsToFile(const char *filename) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
  
  ierr = vsurfaccum.write(filename, NC_FLOAT); CHKERRQ(ierr);
  return 0;
}


