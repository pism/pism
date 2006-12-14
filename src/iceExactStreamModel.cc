// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include <cmath>
#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

#include "iceExactStreamModel.hh"
#include "exactTestI.h"

const PetscScalar IceExactStreamModel::m_schoof = 10; // (pure number)
const PetscScalar IceExactStreamModel::L_schoof = 40e3; // m
const PetscScalar IceExactStreamModel::aspect_schoof = 0.05; // (pure)
const PetscScalar IceExactStreamModel::h0_schoof = aspect_schoof * L_schoof; // = 2000 m
const PetscScalar IceExactStreamModel::B_schoof = 3.7e8; // Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big?
const PetscScalar IceExactStreamModel::p_schoof = 4.0/3.0; // = 1 + 1/n


IceExactStreamModel::IceExactStreamModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {  // do nothing; note derived classes must have constructors
}


void IceExactStreamModel::setflowlawNumber(PetscInt law) {
  flowlawNumber = law;
}


PetscInt IceExactStreamModel::getflowlawNumber() {
  return flowlawNumber;
}


PetscErrorCode IceExactStreamModel::initFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      sometestchosen, testIchosen, inFileSet;
  char            inFile[PETSC_MAX_PATH_LEN], temptestname[20], temp;

  //   "-test I" should already have been chosen, but confirm
  ierr = PetscOptionsGetString(PETSC_NULL, "-test", temptestname, 1, &sometestchosen); CHKERRQ(ierr);
  if (sometestchosen == PETSC_TRUE) {
    temp = temptestname[0];
    if ((temp != 'i') && (temp != 'I')) {
      SETERRQ(1,"IceExactStreamModel only does Test I for now!\n");
    }
  }
  
  /* This switch turns off actual numerical evolution and simply reports the
     exact solution. */
  ierr = PetscOptionsHasName(PETSC_NULL, "-eo", &exactOnly); CHKERRQ(ierr);
  
  /* input file not allowed */
  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    SETERRQ(2,"PISM input file not allowed for initialization of IceExactStreamModel");
  }
  
  ierr = verbPrintf(2,grid.com, 
            "initializing Test I ... \n"); CHKERRQ(ierr);

  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);

  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = grid.rescale(, , ); CHKERRQ(ierr);
  
  ierr = fillinTemps();  CHKERRQ(ierr);
  
  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running EISMINT ROSS ...\n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactStreamModel::fillinTemps() {
  PetscErrorCode      ierr;
  PetscScalar         **Ts, ***T, ***Tb, ***tau;

  // fill in all temps with Ts
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<grid.p->Mz; k++) {
        T[i][j][k] = Ts[i][j];
        tau[i][j][k] = 0.0;
      }
      for (PetscInt k=0; k<grid.p->Mbz; k++)
        Tb[i][j][k] = Ts[i][j];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceExactStreamModel::taucGet(PetscInt i, PetscInt j) const {
  // compare IceCompModel::mapcoords()
  const PetscScalar jfrom0
          = static_cast<PetscScalar>(j) - static_cast<PetscScalar>(grid.p->My - 1)/2.0;
  const PetscScalar y = grid.p->dy * jfrom0;
  const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
  const PetscScalar f = ice.rho * ice.grav * h0_schoof * tan(theta);
  return f * pow(PetscAbs(y / L_schoof), m_schoof);
}


// reimplement basal drag as plastic law [THIS MAY OR MAY NOT WORK!!]
PetscScalar IceExactStreamModel::basalDragx(PetscScalar **u, PetscScalar **v,
                                            PetscInt i, PetscInt j) const {
  return taucGet(i,j) / sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j]));  
}


// ditto
PetscScalar IceExactStreamModel::basalDragy(PetscScalar **u, PetscScalar **v,
                                            PetscInt i, PetscInt j) const {
  return taucGet(i,j) / sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j]));  
}



PetscErrorCode IceExactStreamModel::setBoundaryVels() {
  PetscErrorCode ierr;
  PetscScalar   **ubar, **vbar, **uvbar[2];  // meaning switched!
  
  // set initial velocities in shelf for iteration
  ierr = VecSet(vubar,-200.0 / secpera); CHKERRQ(ierr);
  ierr = VecSet(vvbar,0.0); CHKERRQ(ierr);

  // set boundary condition which will apply to finite difference system:
  //    staggered grid velocities at MASK_SHEET points which neighbor MASK_FLOATING points
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, ubarBC, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vbarBC, &vbar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  for (PetscInt k=0; k < 77; k++) {
    const PetscInt i = kbcGridLoc[0][k];
    const PetscInt j = kbcGridLoc[1][k];
    uvbar[1][i][j-1] = vbar[i][j];
    uvbar[1][i][j] = vbar[i][j];
    uvbar[0][i-1][j] = ubar[i][j];
    uvbar[0][i][j] = ubar[i][j];
  }
  for (PetscInt k=0; k < 22; k++) {
    const PetscInt i = inletGridLoc[0][k];
    const PetscInt j = inletGridLoc[1][k];
    uvbar[1][i][j-1] = vbar[i][j];
    uvbar[1][i][j] = vbar[i][j];
    uvbar[0][i-1][j] = ubar[i][j];
    uvbar[0][i][j] = ubar[i][j];
  }
  ierr = DAVecRestoreArray(grid.da2, ubarBC, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vbarBC, &vbar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactStreamModel::run() {
  PetscErrorCode  ierr;
  PetscInt        pause_time;
  PetscTruth      pause_p, showobsvel, tune;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, &pause_p); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
  "$$$$$      YEAR (+    STEP[$]):     VOL    AREA    MELTF     THICK0     TEMP0\n");
  CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep!

  ierr = setBoundaryVels(); CHKERRQ(ierr);

  // solve model equations 
  useMacayealVelocity = PETSC_TRUE;
  useConstantNuForMacAyeal = PETSC_FALSE;
  useConstantHardnessForMacAyeal = PETSC_TRUE;
  setMacayealEpsilon(0.0);  // don't use this lower bound
//  constantHardnessForMacAyeal = 1.9e8;  // Pa s^{1/3}
  constantHardnessForMacAyeal = 2.22e8;  // Pa s^{1/3}

  regularizingVelocitySchoof = 1.0 / secpera;  // 1 m/a is small velocity for ice stream?
  regularizingLengthSchoof = 1000.0e3;         // (VELOCITY/LENGTH)^2  is very close to 10^-27
  ierr = verbPrintf(5,grid.com,"  [using Schoof regularization constant = %10.5e]\n",
              PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof)); CHKERRQ(ierr);

  ierr = velocityMacayeal(); CHKERRQ(ierr);

  // report on result of computation (i.e. to standard out, to viewers, to Matlab file)
  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  ierr = summary(true,true); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);
  // compare to EXACT ...
  
  if (pause_p == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }

  return 0;
}
