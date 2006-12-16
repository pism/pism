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
#include "exact/exactTestI.h"

const PetscScalar IceExactStreamModel::m_schoof = 10; // (pure number)
const PetscScalar IceExactStreamModel::L_schoof = 40e3; // m
const PetscScalar IceExactStreamModel::aspect_schoof = 0.05; // (pure)
const PetscScalar IceExactStreamModel::H0_schoof = aspect_schoof * L_schoof; // = 2000 m THICKNESS
const PetscScalar IceExactStreamModel::B_schoof = 3.7e8; // Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big?
const PetscScalar IceExactStreamModel::p_schoof = 4.0/3.0; // = 1 + 1/n

// 1 (m/a)^2 is small in basalDrag[x|y] below
const PetscScalar IceExactStreamModel::DEFAULT_PLASTIC_REGULARIZE = 1.0 / secpera; 


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
  PetscTruth      sometestchosen, inFileSet, plasticRegSet;
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
  
  /* This parameter controls regularization of plastic basal sliding law. */
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-plastic_reg", &plastic_regularize,
                              &plasticRegSet); CHKERRQ(ierr);
  if (plasticRegSet == PETSC_FALSE) {
    plastic_regularize = DEFAULT_PLASTIC_REGULARIZE;
  }
  
  /* input file not allowed */
  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    SETERRQ(2,"PISM input file not allowed for initialization of IceExactStreamModel");
  }
  
  ierr = verbPrintf(2,grid.com,"initializing Test I ... \n"); CHKERRQ(ierr);

  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);

  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  
  // set up 8000km by 240km by 3000m grid; note 240km = 6L where L = L_schoof
  // try -Mx 801 -My 25
  ierr = grid.rescale(4000e3, 120e3, 3000); CHKERRQ(ierr);
  
  ierr = fillinTemps();  CHKERRQ(ierr);
  
  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running Test I ...\n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactStreamModel::fillinTemps() {
  PetscErrorCode      ierr;
  const PetscScalar T0 = 263.15;  // completely arbitrary
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);
  ierr = VecSet(vtau, 0.0); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceExactStreamModel::taucGet(PetscInt i, PetscInt j) const {
  // compare IceCompModel::mapcoords()
  const PetscScalar jfrom0
          = static_cast<PetscScalar>(j) - static_cast<PetscScalar>(grid.p->My - 1)/2.0;
  const PetscScalar y = grid.p->dy * jfrom0;
  const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
  const PetscScalar f = ice.rho * ice.grav * H0_schoof * tan(theta);
  return f * pow(PetscAbs(y / L_schoof), m_schoof);
}


// reimplement basal drag as plastic law [THIS MAY OR MAY NOT WORK!!]
PetscScalar IceExactStreamModel::basalDragx(PetscScalar **u, PetscScalar **v,
                                            PetscInt i, PetscInt j) const {
  return taucGet(i,j) / sqrt(plastic_regularize + PetscSqr(u[i][j]) + PetscSqr(v[i][j]));  
}


// ditto
PetscScalar IceExactStreamModel::basalDragy(PetscScalar **u, PetscScalar **v,
                                            PetscInt i, PetscInt j) const {
  return taucGet(i,j) / sqrt(plastic_regularize + PetscSqr(u[i][j]) + PetscSqr(v[i][j]));  
}



PetscErrorCode IceExactStreamModel::setInitStateAndBoundaryVels() {
  PetscErrorCode ierr;
  PetscScalar    Mx = grid.p->Mx, My = grid.p->My;
  PetscScalar    **uvbar[2], **mask, **h, **bed;
  
  // set initial velocities in shelf for iteration
  ierr = VecSet(vubar,0.0); CHKERRQ(ierr);
  ierr = VecSet(vvbar,0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);
  ierr = VecSet(vMask,MASK_DRAGGING); CHKERRQ(ierr);
  ierr = VecSet(vH,H0_schoof); CHKERRQ(ierr);

  // go around edge of grid and set velocity and set SHEET
  // set h, bed everywhere
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk, myu, myv;
      const PetscScalar 
                  ifrom0 = static_cast<PetscScalar>(i)-static_cast<PetscScalar>(Mx - 1)/2.0,
                  jfrom0 = static_cast<PetscScalar>(j)-static_cast<PetscScalar>(My - 1)/2.0;
      const PetscScalar myx = grid.p->dx*ifrom0, myy = grid.p->dy*jfrom0;
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(bed[i][j]), &junk, &myu, &myv); 
      h[i][j] = bed[i][j] + H0_schoof;
      bool edge = (   (i == 0) || (i == 1) || (i == Mx-2) || (i == Mx-1)
                   || (j == 0) || (j == 1) || (j == My-2) || (j == My-1) );
      if (edge) {
        // set boundary condition which will apply to finite difference system:
        // staggered grid velocities at MASK_SHEET points at edges of grid
        mask[i][j] = MASK_SHEET;
        uvbar[0][i-1][j] = myu;
        uvbar[0][i][j] = myu;    // so average onto regular grid point (i,j) has u=myu
        uvbar[1][i][j-1] = myv;
        uvbar[1][i][j] = myv;    // so average onto regular grid point (i,j) has v=myv
      }
    }
  }  
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);    

  return 0;
}


PetscErrorCode IceExactStreamModel::reportErrors() {
//  PetscErrorCode  ierr;
  //NEEDS WRITING!
  return 0;
}

PetscErrorCode IceExactStreamModel::run() {
  PetscErrorCode  ierr;
  PetscInt        pause_time;
  PetscTruth      pause_p;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, &pause_p); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
  "$$$$$      YEAR (+    STEP[$]):     VOL    AREA    MELTF     THICK0     TEMP0\n");
  CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep!

  ierr = fillinTemps(); CHKERRQ(ierr);
  ierr = setInitStateAndBoundaryVels(); CHKERRQ(ierr);

  // solve model equations 
  useMacayealVelocity = PETSC_TRUE;
  useConstantNuForMacAyeal = PETSC_FALSE;
  useConstantHardnessForMacAyeal = PETSC_TRUE;
  setMacayealEpsilon(0.0);  // don't use this lower bound
  constantHardnessForMacAyeal = B_schoof;

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
