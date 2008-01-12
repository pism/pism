// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cstring>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "iceEISModel.hh"


IceEISModel::IceEISModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {  // do nothing; note derived classes must have constructors
  expername = 'A';
  infileused = false;
  flowLawNumber = 0;
}


void IceEISModel::setFlowLawNumber(PetscInt law) {
  flowLawNumber = law;
}


PetscErrorCode IceEISModel::setFromOptions() {
  PetscErrorCode      ierr;

  // apply eismint defaults settings; options may overwrite
  thermalBedrock = PETSC_FALSE;
  useSSAVelocity = PETSC_FALSE;
  isDrySimulation = PETSC_TRUE;
  doGrainSize = PETSC_FALSE;
  enhancementFactor = 1.0;
  includeBMRinContinuity = PETSC_FALSE; // so basal melt does not change 
                                        // computation of vertical velocity

  // optionally allow override of updateHmelt == PETSC_FALSE for EISMINT II
  ierr = PetscOptionsHasName(PETSC_NULL, "-track_Hmelt", &updateHmelt); CHKERRQ(ierr);

  // make bedrock material properties into ice properties
  // (note Mbz=1 is default, but want ice/rock interface segment to 
  // have geothermal flux applied directly to ice)
  bedrock.rho = ice.rho;
  bedrock.k = ice.k;
  bedrock.c_p = ice.c_p;  

  /* This option determines the single character name of EISMINT II experiments:
  "-eisII F", for example.   If not given then do exper A.  */
  char                temp, eisIIexpername[20];
  PetscTruth          EISIIchosen;
  ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", eisIIexpername, 1, &EISIIchosen);
            CHKERRQ(ierr);
  if (EISIIchosen == PETSC_TRUE) {
    temp = eisIIexpername[0];
    if ((temp >= 'a') && (temp <= 'z'))   temp += 'A'-'a';  // capitalize if lower
    if ((temp >= 'A') && (temp <= 'L')) {
      expername = temp;
    } else {
      SETERRQ(1,"option -eisII must have value A, B, C, D, E, F, G, H, I, J, K, or L\n");
    }
  }

  ierr = verbPrintf(2,grid.com, 
              "setting parameters for EISMINT II experiment %c ... \n", 
              expername); CHKERRQ(ierr);
  // EISMINT II specified values for parameters M_max, R_el, T_min, S_b, S_T
  S_b = 1.0e-2 * 1e-3 / secpera;    // Grad of accum rate change
  S_T = 1.67e-2 * 1e-3;           // K/m  Temp gradient
  switch (expername) {
    case 'A':
    case 'E':  // starts from end of A
    case 'G':
    case 'H':
    case 'I':
    case 'K':
      // start with zero ice and:
      M_max = 0.5 / secpera;  // Max accumulation
      R_el = 450.0e3;           // Distance to equil line (accum=0)
      T_min = 238.15;
      break;
    case 'B':
      // supposed to start from end of experiment A and:
      M_max = 0.5 / secpera;
      R_el = 450.0e3;
      T_min = 243.15;
      break;
    case 'C':
    case 'J':
    case 'L':
      // supposed to start from end of experiment A (for C; resp I and K for J and L) and:
      M_max = 0.25 / secpera;
      R_el = 425.0e3;
      T_min = 238.15;
      break;
    case 'D':
      // supposed to start from end of experiment A and:
      M_max = 0.5 / secpera;
      R_el = 425.0e3;
      T_min = 238.15;
      break;
    case 'F':
      // start with zero ice and:
      M_max = 0.5 / secpera;
      R_el = 450.0e3;
      T_min = 223.15;
      break;
    default:
      SETERRQ(999,"\n HOW DID I GET HERE?\n");
  }

  // if user specifies Tmin, Mmax, Sb, ST, Rel, then use that (override above)
  PetscScalar myTmin, myMmax, mySb, myST, myRel;
  PetscTruth  paramSet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Tmin", &myTmin, &paramSet); CHKERRQ(ierr);
  if (paramSet == PETSC_TRUE)     T_min = myTmin;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Mmax", &myMmax, &paramSet); CHKERRQ(ierr);
  if (paramSet == PETSC_TRUE)     M_max = myMmax / secpera;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Sb", &mySb, &paramSet); CHKERRQ(ierr);
  if (paramSet == PETSC_TRUE)     S_b = mySb * 1e-3 / secpera;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ST", &myST, &paramSet); CHKERRQ(ierr);
  if (paramSet == PETSC_TRUE)     S_T = myST * 1e-3;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Rel", &myRel, &paramSet); CHKERRQ(ierr);
  if (paramSet == PETSC_TRUE)     R_el = myRel * 1e3;

  ierr = IceModel::setFromOptions();  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::initFromOptions() {
  PetscErrorCode      ierr;
  PetscTruth          inFileSet, bootFileSet;

  // check if input file was used
  ierr = PetscOptionsHasName(PETSC_NULL, "-if", &inFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-bif", &bootFileSet); CHKERRQ(ierr);
  infileused = ((inFileSet == PETSC_TRUE) || (bootFileSet == PETSC_TRUE));
  
  if (!infileused) { 
    // initialize from EISMINT II formulas
    ierr = verbPrintf(1,grid.com, 
              "initializing EISMINT II experiment %c ... \n", 
              expername); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);

    // following will be part of saved state; not reset if read from file
    // if no inFile then starts with zero ice
    const PetscScalar   G_geothermal   = 0.042;      // J/m^2 s; geo. heat flux
    ierr = VecSet(vh, 0);
    ierr = VecSet(vH, 0);
    ierr = VecSet(vbed, 0);
    ierr = VecSet(vHmelt, 0.0);
    ierr = VecSet(vGhf, G_geothermal);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    ierr = VecSet(vMask, MASK_SHEET);
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);  // no expers use Goldsby-Kohlstedt
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);  // no expers have uplift at start

    // note height of grid must be great enough to handle max thickness
    const PetscScalar   L = 750.0e3;      // Horizontal half-width of grid
    switch (expername) {
      case 'A':
      case 'E':
      case 'I':
        ierr = grid.rescale(L, L, 5000); CHKERRQ(ierr);
        break;
      case 'B':
      case 'C':
      case 'D':
        ierr = grid.rescale(L, L, 4000); CHKERRQ(ierr);
        break;
      case 'F':
        switch (flowLawNumber) {
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
      case 'G':
        ierr = grid.rescale(L, L, 3000); CHKERRQ(ierr);
        break;
      case 'H':
      case 'J':
      case 'K':
      case 'L':
        ierr = grid.rescale(L, L, 5000); CHKERRQ(ierr);
        break;
      default:  
        SETERRQ1(1,"EISMINT II experiment name %c not valid\n",expername);
    }

    if ((expername == 'I') || (expername == 'J')) {
      ierr = generateTroughTopography(); CHKERRQ(ierr);
    } 
    if ((expername == 'K') || (expername == 'L')) {
      ierr = generateMoundTopography(); CHKERRQ(ierr);
    } 

    ierr = initAccumTs(); CHKERRQ(ierr);
    ierr = fillintemps(); CHKERRQ(ierr);

    initialized_p = PETSC_TRUE;
  }

  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);

  if (infileused) {
    ierr = initAccumTs(); CHKERRQ(ierr); // just overwrite accum and Ts with EISMINT II vals
  }
  
  if (infileused && ((expername == 'I') || (expername == 'J'))) { // always regenerate topography
    ierr = generateTroughTopography(); CHKERRQ(ierr); 
  }
  
  ierr = verbPrintf(1,grid.com, "running EISMINT II experiment %c ...\n",expername);
             CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::initAccumTs() {
  PetscErrorCode    ierr;
  PetscScalar       **accum, **Ts;

  // now fill in accum and surface temp
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  PetscScalar cx = grid.p->Lx, cy = grid.p->Ly;
  if (expername == 'E') {  cx += 100.0e3;  cy += 100.0e3;  } // shift center
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // r is distance from center of grid; if E then center is shifted (above)
      const PetscScalar r = sqrt( PetscSqr(-cx + grid.p->dx*i) + PetscSqr(-cy + grid.p->dy*j) );
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
  PetscScalar         **Ts;
  const PetscScalar   G_geothermal   = 0.042; // J/m^2 s; geo. heat flux; only matters if
                                              // (nonstandard case) bedrock is present

  // fill in all ice temps with Ts and then have bedrock (if present despite EISMINT
  //   standard) temperatures reflect default geothermal rate
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = T3.needAccessToVals(); CHKERRQ(ierr);
  ierr = Tb3.needAccessToVals(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.setToConstantColumn(i,j,Ts[i][j]); CHKERRQ(ierr);
      ierr = bootstrapSetBedrockColumnTemp(i,j,Ts[i][j],G_geothermal); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = T3.needAccessToVals(); CHKERRQ(ierr);
  ierr = Tb3.needAccessToVals(); CHKERRQ(ierr);

  // communicate T because it will be horizontally differentiated
  ierr = T3.beginGhostComm(); CHKERRQ(ierr);
  ierr = T3.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::generateTroughTopography() {
  PetscErrorCode  ierr;
  // computation based on
  //    http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f
  // by Tony Payne, 6 March 1997
  
  const PetscScalar    b0 = 1000.0;  // plateau elevation
  const PetscScalar    L = 750.0e3;  // half-width of computational domain
  const PetscScalar    w = 200.0e3;  // trough width
  const PetscScalar    slope = b0/L;
  const PetscScalar    dx = grid.p->dx, dy = grid.p->dy;
  const PetscScalar    dx61 = (2*L) / 60; // = 25.0e3
  PetscScalar          topg, **b;

  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar nsd = i * dx, ewd = j *dy;
      if (    (nsd >= (27 - 1) * dx61) && (nsd <= (35 - 1) * dx61)
           && (ewd >= (31 - 1) * dx61) && (ewd <= (61 - 1) * dx61) ) {
        topg = 1000.0 - PetscMax(0.0, slope * (ewd - L) * cos(pi * (nsd - L) / w));
      } else {
        topg = 1000.0;
      }
      b[i][j] = topg;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  // communicate b because it will be horizontally differentiated
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
               "trough bed topography stored by IceEISModel::generateTroughTopography()\n");
               CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::generateMoundTopography() {
  PetscErrorCode  ierr;
  // computation based on
  //    http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f
  // by Tony Payne, 6 March 1997
  
  const PetscScalar    slope = 250.0;
  const PetscScalar    w = 150.0e3;  // mound width
  const PetscScalar    dx = grid.p->dx, dy = grid.p->dy;
  PetscScalar          topg, **b;

  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar nsd = i * dx, ewd = j *dy;
      topg = PetscAbs(slope * sin(pi * ewd / w) + slope * cos(pi * nsd / w));
      b[i][j] = topg;
    }
  }
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  // communicate b because it will be horizontally differentiated
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "mound bed topography stored by IceEISModel::generateTroughTopography()\n");
           CHKERRQ(ierr);
  return 0;
}


// reimplement IceModel::basalVelocity() which is virtual; basalVelocity() is 
// for SIA regions (MASK_SHEET), and it is called within IceModel::velocitySIAStaggered)
PetscScalar IceEISModel::basalVelocity(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {

  const PetscScalar  Bfactor = 1e-3 / secpera; // units m s^-1 Pa^-1
  const PetscScalar  eismintII_temp_sliding = 273.15;
  
  if (expername == 'G') {
      return Bfactor * ice.rho * grav * H; 
  } else if (expername == 'H') {
      if (T + ice.beta_CC_grad * H > eismintII_temp_sliding) {
        return Bfactor * ice.rho * grav * H; // ditto case G
      } else {
        return 0.0;
      }
  }  
  return 0.0;  // zero sliding for other tests
}

