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
#include "pccoupler.hh"


/******************* VIRTUAL BASE CLASS:  PISMClimateCoupler ********************/

PISMClimateCoupler::PISMClimateCoupler() {
}


PISMClimateCoupler::~PISMClimateCoupler() { 
}


//! User like IceModel will call this to set the grid, so that PISMClimateCoupler can use it.
PetscErrorCode PISMClimateCoupler::setGrid(IceGrid* g) {
  grid = g;
  return 0;
}


PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(const char *filename) {
  SETERRQ(1,"not implemented");
  return 0;
}


/******************* ATMOSPHERE:  PISMAtmosphereCoupler ********************/

/*!
We set attributes here but do not allocate space.  A derived class must
allocate.  (Attempting to get an array pointer to vsurfmassflux, vsurftemp
will fail in this base class.)
 */
PISMAtmosphereCoupler::PISMAtmosphereCoupler() : PISMClimateCoupler() {

  // check whether standard_name is on CF table
  //   http://cf-pcmdi.llnl.gov/documents/cf-standard-names/current/cf-standard-name-table.html
  // i.e. COMPARE TO CURRENT VERSION IN IceModel: vAccum, vTs
  
  // arguments to set_attrs() are: pism_intent, long_name, units, standard_name;

  vsurfmassflux.set_attrs(
           "climatecoupler",
           "ice equivalent net mass balance at surface of ice",
	   "m/s",
	   NULL); 	
  vsurfmassflux.set_glaciological_units("m year-1", secpera);
  vsurftemp.set_attrs(
           "climatecoupler",
           "absolute ice temperature at surface of ice",
	   "K",
	   NULL);
}


PISMAtmosphereCoupler::~PISMAtmosphereCoupler() {
  vsurfmassflux.destroy();
  vsurftemp.destroy();
}


PetscErrorCode PISMAtmosphereCoupler::writeCouplingFieldsToFile(const char *filename) {
  SETERRQ(1,"not implemented");
  return 0;
}


/*******************  OCEAN:  PISMOceanCoupler ********************/

PISMOceanCoupler::PISMOceanCoupler() : PISMClimateCoupler() {
}

