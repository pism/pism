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


PISMClimateCoupler::PISMClimateCoupler() {
}


PISMClimateCoupler::~PISMClimateCoupler() { 
}


PetscErrorCode PISMClimateCoupler::setGrid(IceGrid* g) {
  grid = g;
  return 0;
}


PetscErrorCode PISMClimateCoupler::init() {
  SETERRQ(1,"not implemented");
  return 0;
}

PetscErrorCode PISMClimateCoupler::writeCouplingFieldsToFile(const char *filename) {
  SETERRQ(1,"not implemented");
  return 0;
}


PISMAtmosphereCoupler::PISMAtmosphereCoupler() : PISMClimateCoupler() {
}


PISMOceanCoupler::PISMOceanCoupler() : PISMClimateCoupler() {
}

