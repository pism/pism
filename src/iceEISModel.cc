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
#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

#include "iceEISModel.hh"


IceEISModel::IceEISModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {  // do nothing; note derived classes must have constructors
}


void IceEISModel::setExperName(char name) {
  expername = name;
}


char IceEISModel::getExperName() {
  return expername;
}


void IceEISModel::setflowlawNumber(PetscInt law) {
  flowlawNumber = law;
}


PetscInt IceEISModel::getflowlawNumber() {
  return flowlawNumber;
}


PetscErrorCode IceEISModel::setExperNameFromOptions() {
  PetscErrorCode      ierr;
  char                temp, eisIIexpername[20];
  PetscTruth          EISIIchosen;

  /* note EISMINT I is NOT worth implementing; for fixed margin isothermal 
     tests do "pismv -test A" or "pismv -test E"; for moving margin isothermal
     tests do "pismv -test B" or "-test C" or "-test D" or "-test H" */

  /* This option determines the single character name of EISMINT II experiments:
  "-eisII F", for example. */
  ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", eisIIexpername, 1, &EISIIchosen); CHKERRQ(ierr);

//   This option chooses EISMINT ROSS, i.e. from the paper
//     MacAyeal and five others (1996). "An ice-shelf model test based on the Ross ice shelf,"
//     Ann. Glaciol. 23, 46-51
//  ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &ROSSchosen); CHKERRQ(ierr);

  if (EISIIchosen == PETSC_TRUE) {
    temp = eisIIexpername[0];
    if ((temp >= 'a') && (temp <= 'z'))   temp += 'A'-'a';  // capitalize if lower
    setExperName(temp);
//  } else if (ROSSchosen == PETSC_TRUE) {
//    setExperName('3');
//    SETERRQ(1, "EISMINT-ROSS NOT IMPLEMENTED!");
  } else { // set a default: EISMINT II experiment A
    setExperName('A');
  }

  return 0;
}


PetscErrorCode IceEISModel::initFromOptions() {
  PetscErrorCode      ierr;
  char                inFile[PETSC_MAX_PATH_LEN];
  const PetscScalar   G_geothermal   = 0.042;      // J/m^2 s; geo. heat flux
  const PetscScalar   L              = 750e3;      // Horizontal extent of grid

  ierr = setExperNameFromOptions(); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, 
              "initializing EISMINT II experiment %c ... \n", 
              getExperName()); CHKERRQ(ierr);
    ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);
    
    // following will be part of saved state; not reset if read from file
    // if no inFile then starts with zero ice
    ierr = VecSet(vh, 0);
    ierr = VecSet(vH, 0);
    ierr = VecSet(vbed, 0);
    ierr = VecSet(vHmelt, 0.0);
    ierr = VecSet(vGhf, G_geothermal);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    ierr = VecSet(vMask, MASK_SHEET);
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);  // no expers use Goldsby-Kohlstedt
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);  // no expers have uplift at start

    // next block checks if experiment is implemented
    // note height of grid must be great enough to handle max thickness
    switch (getExperName()) {
      case 'A':
        ierr = grid.rescale(L, L, 4500); CHKERRQ(ierr);
        break;
      case 'B':
      case 'C':
      case 'D':
        ierr = grid.rescale(L, L, 4000); CHKERRQ(ierr);
        break;
      case 'F':
        switch (getflowlawNumber()) {
          case 0:
          case 3:
            ierr = grid.rescale(L, L, 5000); CHKERRQ(ierr);
            break;
          case 1:
          case 4:
          case 5:
            ierr = grid.rescale(L, L, 6000); CHKERRQ(ierr);
            break;
          case 2:
            ierr = grid.rescale(L, L, 7000); CHKERRQ(ierr);
            break;
          default:  SETERRQ(1,"should not reach here (switch for rescale)\n");
        }
        break;
      default:  
        SETERRQ(1,"ERROR: desired EISMINT II experiment NOT IMPLEMENTED\n");
    }
  }

  ierr = applyDefaultsForExperiment(); CHKERRQ(ierr);
  ierr = initAccumTs(); CHKERRQ(ierr);
  if (inFileSet == PETSC_FALSE) {
    ierr = fillintemps(); CHKERRQ(ierr);
  }
  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = PetscPrintf(grid.com, "running EISMINT II experiment %c ...\n",getExperName());
             CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEISModel::applyDefaultsForExperiment() {

  setThermalBedrock(PETSC_FALSE);
  setUseMacayealVelocity(PETSC_FALSE); // NOT VALID FOR EISMINT-ROSS!!
  setIsDrySimulation(PETSC_TRUE); // NOT VALID FOR EISMINT-ROSS!!
  setDoGrainSize(PETSC_FALSE);
  setEnhancementFactor(1.0);
  setIncludeBMRinContinuity(PETSC_FALSE);

  // make bedrock material properties into ice properties
  // (note Mbz=0 is default, but want ice/rock interface segment to see all ice)
  bedrock.rho = ice.rho;
  bedrock.k = ice.k;
  bedrock.c_p = ice.c_p;  

  return 0;
}


PetscErrorCode IceEISModel::initAccumTs() {
  PetscScalar       M_max,R_el,T_min;
  PetscErrorCode    ierr;
  PetscScalar       **accum, **Ts;

  // EISMINT II values:
  const PetscScalar S_b = 1e-2 * 1e-3 / secpera;    // Grad of accum rate change
  PetscScalar       S_T = 1.67e-2 * 1e-3;           // K/m  Temp gradient

  switch (getExperName()) {
    case 'A':
      // start with zero ice and:
      M_max = 0.5 / secpera;  // Max accumulation
      R_el = 450e3;           // Distance to equil line (accum=0)
      T_min = 238.15;
      break;
    case 'B':
      // supposed to start from end of experiment A and:
      M_max = 0.5 / secpera;
      R_el = 450e3;
      T_min = 243.15;
      break;
    case 'C':
      // supposed to start from end of experiment A and:
      M_max = 0.25 / secpera;
      R_el = 425e3;
      T_min = 238.15;
      break;
    case 'D':
      // supposed to start from end of experiment A and:
      M_max = 0.5 / secpera;
      R_el = 425e3;
      T_min = 238.15;
      break;
    case 'F':
      // start with zero ice and:
      M_max = 0.5 / secpera;
      R_el = 450e3;
      T_min = 223.15;
      break;
    default:
      SETERRQ(1,"\n experiment name unknown in IceEISModel::initAccumTs()\n");
  }

  // now fill in accum and surface temp
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // r is distance from center of grid
      const PetscScalar r = sqrt( PetscSqr(-grid.p->Lx + grid.p->dx*i)
                                  + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
      // set accumulation
      accum[i][j] = PetscMin(M_max, S_b * (R_el-r));  // formula (7) in (Payne et al 2000)
      // set surface temperature
      Ts[i][j] = T_min + S_T * r;                 // formula (8) in (Payne et al 2000)
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::fillintemps() {
  PetscErrorCode      ierr;
  PetscScalar         **Ts, ***T, ***Tb;

  // fill in all temps with Ts
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<grid.p->Mz; k++)
        T[i][j][k] = Ts[i][j];
      for (PetscInt k=0; k<grid.p->Mbz; k++)
        Tb[i][j][k] = Ts[i][j];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  return 0;
}


// reimplement IceModel::basal()
PetscScalar IceEISModel::basal(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {
  // note this version ignors mu
  
  return 0.0;  // zero sliding for other tests; EISMINT II exper G,H will need sliding
}

