// Copyright (C) 2009-2010 Andy Aschwanden and Ed Bueler
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

#include "../../base/grid.hh"
#include "../../base/materials.hh"
#include "../../base/iceModel.hh"
#include "iceScandModel.hh"


IceScandModel::IceScandModel(IceGrid &g) : IceEISModel(g) {
  expername = 'S';
}


PetscErrorCode IceScandModel::set_expername_from_options() {
  PetscErrorCode ierr;
  // do nothing except stop if -eisII is set; default expername already set to S
  PetscTruth eisIISet;
  ierr = check_option("-eisII", eisIISet); CHKERRQ(ierr);
  if (eisIISet == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
      "IceScandModel ERROR: -eisII option not used because does '-eisII S' always\n");
      CHKERRQ(ierr);
    PetscEnd();
  }
  return 0;
}


PetscErrorCode IceScandModel::setFromOptions() {
  PetscErrorCode      ierr;

  ierr = IceEISModel::setFromOptions(); CHKERRQ(ierr);
  
  // start with zero ice and:
  M_max = 0.5 / secpera;  // Max accumulation
  R_el = 450.0e3;           // Distance to equil line (accum=0)
  T_min = 238.15;
  T_max = 273.15;

  R_cts = 100.0e3; // position where transition from temperate to cold
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Rcts", &R_cts, NULL); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceScandModel::init_couplers() {
  PetscErrorCode      ierr;

  // kill off sub-ice shelf initialization blurb; we don't have ocean interaction
  if (oceanPCC != PETSC_NULL) {
    oceanPCC->reportInitializationToStdOut = false;
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
  }

  ierr = IceModel::init_couplers(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
    "  setting surface mass balance and surface temperature variables for IceScandModel ...\n");
    CHKERRQ(ierr);

  PetscTruth i_set;
  char filename[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-i",
			       filename, PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
  if (i_set) {
    ierr = verbPrintf(2,grid.com,
      "  (values from file %s ignored)\n", filename); CHKERRQ(ierr);
  }

  IceModelVec2 *pccsmf, *pccTs;

  if (atmosPCC != PETSC_NULL) {
    // call sets pccsmf to point to IceModelVec2 with current surface massflux
    ierr = atmosPCC->updateSurfMassFluxAndProvide(
              grid.year, 0.0, pccsmf); CHKERRQ(ierr);
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, 0.0, pccTs); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }

  ierr =  pccTs->begin_access(); CHKERRQ(ierr);
  ierr = pccsmf->begin_access(); CHKERRQ(ierr);
  PetscScalar cx = grid.Lx, cy = grid.Ly;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // r is distance from center of grid; if E then center is shifted (above)
      const PetscScalar r = sqrt( PetscSqr(-cx + grid.dx*i)
                                  + PetscSqr(-cy + grid.dy*j) );
      // set accumulation from formula (7) in (Payne et al 2000)
      (*pccsmf)(i,j) = PetscMin(M_max, S_b * (R_el-r));
      // set surface temperature
      // simplest possible Scandinavian-type upper surface boundary condition
      // could be replace with more elaborate formula
      if (r <= R_cts)  (*pccTs)(i,j) = T_max;
      else             (*pccTs)(i,j) = T_min;
    }
  }
  ierr = pccsmf->end_access(); CHKERRQ(ierr);
  ierr =  pccTs->end_access(); CHKERRQ(ierr);
  return 0;
}

