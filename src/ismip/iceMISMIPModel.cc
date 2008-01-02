// Copyright (C) 2008 Ed Bueler
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
#include "iceMISMIPModel.hh"



class MISMIPBasalType : public ViscousBasalType {
public:
  MISMIPBasalType(const PetscScalar m, const PetscScalar C);
  virtual PetscErrorCode printInfo(const int thresh, MPI_Comm com);
  virtual PetscScalar    drag(PetscScalar coeff, PetscScalar tauc, PetscScalar vx, PetscScalar vy);
  PetscScalar m_MISMIP, C_MISMIP;
};


MISMIPBasalType::MISMIPBasalType(const PetscScalar m, const PetscScalar C) {
  m_MISMIP = m;
  C_MISMIP = C;
}


PetscErrorCode MISMIPBasalType::printInfo(const int thresh, MPI_Comm com) {
  PetscErrorCode ierr;
  ierr = verbPrintf(thresh, com, 
    "Using  tau_b = - C |u|^{m-1} u  model for MISMIP, with m=%5.4f and C=%5.4f.\n",
    m_MISMIP, C_MISMIP); CHKERRQ(ierr);
  return 0;
}


PetscScalar MISMIPBasalType::drag(PetscScalar coeff, PetscScalar tauc,
                                   PetscScalar vx, PetscScalar vy) {
  const PetscScalar magsliding = sqrt(PetscSqr(vx) + PetscSqr(vy));
  return C_MISMIP * pow(magsliding, m_MISMIP - 1.0);
}



IceMISMIPModel::IceMISMIPModel(IceGrid &g, IceType &i) : IceModel(g, i) {
  // only call parent's constructor; do all classes need constructors?
}


PetscErrorCode IceMISMIPModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  
  return 0;
}


PetscErrorCode IceMISMIPModel::initFromOptions() {
  PetscErrorCode ierr;
  char inFile[PETSC_MAX_PATH_LEN];
  PetscTruth inFileSet, bootFileSet;
  
  ierr = IceModel::initFromOptions(PETSC_FALSE); CHKERRQ(ierr);  // wait on init hook; possible regridding!

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootFileSet); CHKERRQ(ierr);
  
  if (inFileSet == PETSC_TRUE) {
    if (bootFileSet) {
      ierr = verbPrintf(1, grid.com, "WARNING: -bif and -if given; using -if\n"); CHKERRQ(ierr);
    }
  } 
  
  if (!isInitialized()) {
    SETERRQ(1, "IceMISMIPModel has not been initialized.\n");
  }

  ierr = afterInitHook(); CHKERRQ(ierr);  // note regridding can happen here
  return 0;
}


PetscErrorCode IceMISMIPModel::additionalAtStartTimestep() {
  // PetscErrorCode ierr;

  return 0;
}




