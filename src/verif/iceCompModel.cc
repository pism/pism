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

#include <cmath>
#include <cstring>
#include <petscda.h>

#include "exactTestsABCDE.h"
#include "exactTestsFG.h" 
#include "exactTestH.h" 
#include "exactTestL.h" 
#include "../num/extrasGSL.hh"

#include "iceCompModel.hh"

PetscScalar IceCompModel::ablationRateOutside = 0.02; // m/a

IceCompModel::IceCompModel(IceGrid &g, ThermoGlenArrIce &i, const char mytest)
  : IceModel(g, i), tgaIce(i) {
  
  // note lots of defaults are set by the IceModel constructor
  
  testname = mytest;
  
  // Override some defaults from parent class
  enhancementFactor = 1.0;
  f = tgaIce.rho / bedrock.rho;  // for simple isostasy
  
  // defaults for verification
  exactOnly = PETSC_FALSE;
  compVecsCreated = PETSC_FALSE;
  compViewersCreated = PETSC_FALSE;
  vHexactLCreated = PETSC_FALSE;

  // set values of flags in run() 
  doMassConserve = PETSC_TRUE;
  useSSAVelocity = PETSC_FALSE;
  includeBMRinContinuity = PETSC_FALSE;
  doPlasticTill = PETSC_FALSE;
  doPDD = PETSC_FALSE;
  doGrainSize = PETSC_FALSE;
  if (testname == 'H') {
    doBedDef = PETSC_TRUE;
    doBedIso = PETSC_TRUE;
  } else
    doBedDef = PETSC_FALSE;  
  if ((testname == 'F') || (testname == 'G') || (testname == 'K')) {
    doTemp = PETSC_TRUE;
    globalMinAllowedTemp = 0.0;  // essentially turn off run-time reporting of extremely
    maxLowTempCount = 1000000;   // low computed temperatures; *they will be reported
                                 // as errors* anyway
  } else
    doTemp = PETSC_FALSE; 
  if ((testname == 'A') || (testname == 'E')) {
    isDrySimulation = PETSC_TRUE;
    doOceanKill = PETSC_TRUE;
  } else {
    isDrySimulation = PETSC_TRUE;
    doOceanKill = PETSC_FALSE;
  }

  // special considerations for K wrt thermal bedrock and pressure-melting
  bedrock_is_ice_forK = PETSC_FALSE;
  if (testname == 'K') {
    thermalBedrock = PETSC_TRUE;
    allowAboveMelting = PETSC_FALSE;
    reportHomolTemps = PETSC_TRUE;
  } else {
    thermalBedrock = PETSC_FALSE;
    // note temps are generally allowed to go above pressure melting in verify
    allowAboveMelting = PETSC_TRUE;
    reportHomolTemps = PETSC_FALSE;
    // now make bedrock have same material properties as ice
    // (note Mbz=1 also, by default, but want ice/rock interface to see
    // pure ice from the point of view of applying geothermal boundary
    // condition, especially in tests F and G)
    bedrock.rho = tgaIce.rho;
    bedrock.c_p = tgaIce.c_p;
    bedrock.k = tgaIce.k;
  }
}


IceCompModel::~IceCompModel() {
  if (compViewersCreated == PETSC_TRUE) {
    destroyCompViewers();
    compViewersCreated = PETSC_FALSE;
  }
  if (compVecsCreated == PETSC_TRUE) {
    destroyCompVecs();
    compVecsCreated = PETSC_FALSE;
  }
  if (vHexactLCreated == PETSC_TRUE) {
    VecDestroy(vSigmaComp);
    vHexactLCreated = PETSC_FALSE;
  }
}


PetscErrorCode IceCompModel::initFromOptions() {
  //  do this only after IceCompModel::setFromOptions()
  PetscErrorCode ierr;
  PetscTruth inFileSet;
  char inFile[PETSC_MAX_PATH_LEN];

  if ( (testname == 'H') && ((doBedDef == PETSC_FALSE) || (doBedIso == PETSC_FALSE)) ) {
    ierr = verbPrintf(1,grid.com, "[IceCompModel WARNING: Test H should be run with option  ");
    ierr = verbPrintf(1,grid.com, "-bed_def_iso  for the reported errors to be correct.]\n");
    CHKERRQ(ierr);
  }

  /* This switch turns off actual numerical evolution and simply reports the
     exact solution. */
  ierr = PetscOptionsHasName(PETSC_NULL, "-eo", &exactOnly); CHKERRQ(ierr);

  /* This switch turns changes Test K to make material properties for bedrock
     the same as for the ice. */
  PetscTruth biiSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-bedrock_is_ice", &biiSet); CHKERRQ(ierr);
  if (biiSet == PETSC_TRUE) {
    if (testname == 'K') {
      ierr = verbPrintf(1,grid.com,
              "setting material properties of bedrock to those of ice in Test K\n");
              CHKERRQ(ierr);
      bedrock.rho = tgaIce.rho;
      bedrock.c_p = tgaIce.c_p;
      bedrock.k = tgaIce.k;
      bedrock_is_ice_forK = PETSC_TRUE;
    } else {
      ierr = verbPrintf(1,grid.com,
              "WARNING: option -bedrock_is_ice ignored; only applies to Test K\n");
              CHKERRQ(ierr);
    }
  }

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);

  if (inFileSet == PETSC_FALSE) {
    ierr = verbPrintf(2,grid.com, "initializing Test %c ...\n",testname);  CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);
    switch (testname) {
      case 'A':
      case 'E':
        // use 1600km by 1600km by 4000m rectangular domain
        ierr = grid.rescale(800e3, 800e3, 4000); CHKERRQ(ierr);
        break;
      case 'B':
        // use 2400km by 2400km by 4000m rectangular domain
        ierr = grid.rescale(1200e3, 1200e3, 4000); CHKERRQ(ierr);
        break;
      case 'C':
      case 'D':
        // use 2000km by 2000km by 4000m rectangular domain
        ierr = grid.rescale(1000e3, 1000e3, 4000); CHKERRQ(ierr);
        break;
      case 'F':
      case 'G':
      case 'L':
        // use 1800km by 1800km by 4000m rectangular domain
        ierr = grid.rescale(900e3, 900e3, 4000); CHKERRQ(ierr);
        break;
      case 'H':
        // use 1500km by 1500km by 4000m rectangular domain
        ierr = grid.rescale(1500e3, 1500e3, 4000); CHKERRQ(ierr);
        break;
      case 'K':
        // use 2000km by 2000km by 4000m rectangular domain, but make truely periodic
        ierr = grid.rescale(1000e3, 1000e3, 4000, PETSC_TRUE); CHKERRQ(ierr);
        break;
      default:  SETERRQ(1,"IceCompModel ERROR : desired test not implemented\n");
    }

    // none use Goldsby-Kohlstedt or do age calc
    ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS); CHKERRQ(ierr);
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    // all have no uplift or Hmelt
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);
    ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);
    ierr = VecSet(vbasalMeltRate, 0.0); CHKERRQ(ierr);

    ierr = createCompVecs(); CHKERRQ(ierr);

    if (yearsStartRunEndDetermined == PETSC_FALSE) {
      ierr = setStartRunEndYearsFromOptions(PETSC_FALSE);  CHKERRQ(ierr);
    }

    switch (testname) {
      case 'A':
      case 'B':
      case 'C':
      case 'D':
      case 'E':
      case 'H':
        ierr = initTestABCDEH(); CHKERRQ(ierr);
        break;
      case 'F':
      case 'G':
        ierr = initTestFG(); CHKERRQ(ierr);  // see iCMthermo.cc
        break;
      case 'K':
        ierr = initTestK(); CHKERRQ(ierr);  // see iCMthermo.cc
        break;
      case 'L':
        ierr = initTestL(); CHKERRQ(ierr);
        break;
      default:  SETERRQ(1,"Desired test not implemented by IceCompModel.\n");
    }

    initialized_p = PETSC_TRUE;
  }

  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);

  if (inFileSet == PETSC_TRUE) {
    ierr = createCompVecs(); CHKERRQ(ierr);
  }

  ierr = createCompViewers(); CHKERRQ(ierr);

  if (exactOnly == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com, "!!EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION!!\n");
             CHKERRQ(ierr);
  }
  return 0;
}


void IceCompModel::mapcoords(const PetscInt i, const PetscInt j,
                             PetscScalar &x, PetscScalar &y, PetscScalar &r) {
  // compute x,y,r on grid from i,j
  PetscScalar ifrom0, jfrom0;

  ifrom0=static_cast<PetscScalar>(i)-static_cast<PetscScalar>(grid.p->Mx - 1)/2.0;
  jfrom0=static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.p->My - 1)/2.0;
  x=grid.p->dx*ifrom0;
  y=grid.p->dy*jfrom0;
  r = sqrt(PetscSqr(x) + PetscSqr(y));
}


// reimplement IceModel::basalVelocity()
PetscScalar IceCompModel::basalVelocity(const PetscScalar xIN, const PetscScalar yIN,
                                        const PetscScalar H, const PetscScalar T,
                                        const PetscScalar alpha, const PetscScalar muIN) {
  // note: ignors T and muIN

  if (testname == 'E') {
    //PetscErrorCode  ierr = PetscPrintf(grid.com, 
    //        "   [IceCompModel::basal called with:   x=%f, y=%f, H=%f, T=%f, alpha=%f]\n",
    //        xIN,yIN,H,alpha);  CHKERRQ(ierr);
    const PetscScalar pi = 3.14159265358979;
    const PetscScalar r1 = 200e3, r2 = 700e3,   /* define region of sliding */
                      theta1 = 10 * (pi/180), theta2 = 40 * (pi/180);
    const PetscScalar x = fabs(xIN), y = fabs(yIN);
    const PetscScalar r = sqrt(x * x + y * y);
    PetscScalar       theta;
    if (x < 1.0)
      theta = pi / 2.0;
    else
      theta = atan(y / x);
  
    if ((r > r1) && (r < r2) && (theta > theta1) && (theta < theta2)) {
      // now INSIDE sliding region
      const PetscScalar rbot = (r2 - r1) * (r2 - r1),
                        thetabot = (theta2 - theta1) * (theta2 - theta1);
      const PetscScalar mu_max = 2.5e-11; /* Pa^-1 m s^-1; max sliding coeff */
      PetscScalar muE = mu_max * (4.0 * (r - r1) * (r2 - r) / rbot) 
                               * (4.0 * (theta - theta1) * (theta2 - theta) / thetabot);
      return muE * tgaIce.rho * grav * H;
    } else
      return 0.0;
  } else
    return 0.0;  // zero sliding for other tests
}


PetscErrorCode IceCompModel::initTestABCDEH() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, **H, **accum, **mask, dummy1, dummy2, dummy3;
  const PetscScalar LforAE = 750e3; // m

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = -tgaIce.Q() / (gasConst_R * log(A0/tgaIce.A()));
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);
  
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  muSliding = 0.0;  // note reimplementation of basalVelocity()

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  if ((testname == 'A') || (testname == 'E')) {
    ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  }
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'A':
          exactA(r,&H[i][j],&accum[i][j]);
          if (r >= LforAE)
            mask[i][j] = MASK_FLOATING_OCEAN0;
          break;
        case 'B':
          exactB(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'C':
          exactC(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'D':
          exactD(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'E':
          exactE(xx,yy,&H[i][j],&accum[i][j],&dummy1,&dummy2,&dummy3);
          if (r >= LforAE)
            mask[i][j] = MASK_FLOATING_OCEAN0;
          break;
        case 'H':
          exactH(f,grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        default:  SETERRQ(1,"test must be A, B, C, D, E, or H");
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  if ((testname == 'A') || (testname == 'E')) {
    ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  }

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);

  if (testname == 'H') {
    ierr = VecCopy(vH,vh); CHKERRQ(ierr);
    ierr = VecScale(vh,1-f); CHKERRQ(ierr);
    ierr = VecCopy(vH,vbed); CHKERRQ(ierr);
    ierr = VecScale(vbed,-f); CHKERRQ(ierr);
  } else {  // flat bed case otherwise
    ierr = VecCopy(vH,vh); CHKERRQ(ierr);
    ierr = VecSet(vbed, 0.0); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceCompModel::initTestL() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, **H, **accum, **bed;

  if (testname != 'L')  { SETERRQ(1,"test must be 'L'"); }
  
  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = 1.0e-16/secpera;    // = 3.17e-24  1/(Pa^3 s);  (EISMINT value) flow law parameter
  T0 = -tgaIce.Q() / (gasConst_R * log(A0/tgaIce.A()));
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);
  
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  muSliding = 0.0;  // note reimplementation of basalVelocity()

  // setup to evaluate test L; requires solving an ODE numerically
  const PetscInt  MM = grid.xm * grid.ym;
  double          *rr, *HH, *bb, *aa;
  int             *ia, *ja;
  rr = new double[MM];  
  HH = new double[MM];  bb = new double[MM];  aa = new double[MM];
  ia = new int[MM];  ja = new int[MM];

  for (PetscInt i = 0; i < grid.xm; i++) {
    for (PetscInt j = 0; j < grid.ym; j++) {
      const PetscInt  k = i * grid.ym + j;
      PetscScalar  junkx, junky;
      mapcoords(i + grid.xs, j + grid.ys, junkx, junky, rr[k]);
      rr[k] = - rr[k];
      ia[k] = i + grid.xs;  ja[k] = j + grid.ys;
    }
  }

  heapsort_double_2indfollow(rr,ia,ja,MM);  // sorts into ascending;  O(MM log MM) time
  for (PetscInt k = 0; k < MM; k++)   rr[k] = -rr[k];   // now descending

  // get soln to test L at these points; solves ODE only once (on each processor)
  ierr = exactL_list(rr, MM, HH, bb, aa);  CHKERRQ(ierr);
  
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  for (PetscInt k = 0; k < MM; k++) {
    const PetscInt i = ia[k],  j = ja[k];
    H[i][j] = HH[k];
    bed[i][j] = bb[k];
    accum[i][j] = aa[k];
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);

  delete [] rr;  delete [] HH;  delete [] bb;  delete [] aa;  delete [] ia;  delete [] ja; 

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);

  // set surface to H+b
  ierr = VecWAXPY(vh,1.0,vH,vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);

  // store copy of vH for "-eo" runs and for evaluating geometry errors
  ierr = VecDuplicate(vh, &vHexactL); CHKERRQ(ierr);
  ierr = VecCopy(vH, vHexactL); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::getCompSourcesTestCDH() {
  PetscErrorCode  ierr;
  PetscScalar     **accum, dummy;

  // before flow step, set accumulation from exact values;
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'C':
          exactC(grid.p->year*secpera,r,&dummy,&accum[i][j]);
          break;
        case 'D':
          exactD(grid.p->year*secpera,r,&dummy,&accum[i][j]);
          break;
        case 'H':
          exactH(f,grid.p->year*secpera,r,&dummy,&accum[i][j]);
          break;
        default:  SETERRQ(1,"testname must be C, D, or H");
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestABCDH() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **accum;

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'A':
          exactA(r,&H[i][j],&accum[i][j]);
          break;
        case 'B':
          exactB(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'C':
          exactC(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'D':
          exactD(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        case 'H':
          exactH(f,grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
          break;
        default:  
          SETERRQ(1,"test must be A, B, C, D, or H");
          break;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);

  if (testname == 'H') {
    ierr = VecCopy(vH,vh); CHKERRQ(ierr);
    ierr = VecScale(vh,1-f); CHKERRQ(ierr);
    ierr = VecCopy(vH,vbed); CHKERRQ(ierr);
    ierr = VecScale(vbed,-f); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(vH,vh); CHKERRQ(ierr);
  }
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestE() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **accum, **ub, **vb, dummy;

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      exactE(xx,yy,&H[i][j],&accum[i][j],&dummy,&ub[i][j],&vb[i][j]);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = VecCopy(vH,vh); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vub, INSERT_VALUES, vub); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vub, INSERT_VALUES, vub); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vvb, INSERT_VALUES, vvb); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vvb, INSERT_VALUES, vvb); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestL() {
  PetscErrorCode  ierr;

  ierr = DALocalToLocalBegin(grid.da2, vHexactL, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vHexactL, INSERT_VALUES, vH); CHKERRQ(ierr);

  ierr = VecWAXPY(vh,1.0,vH,vbed); CHKERRQ(ierr);  // h = H + bed = 1 * H + bed
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);

  // note bed was filled at initialization and hasn't changed
  return 0;
}


PetscErrorCode IceCompModel::computeGeometryErrors(
      PetscScalar &gvolexact, PetscScalar &gareaexact, PetscScalar &gdomeHexact,
      PetscScalar &volerr, PetscScalar &areaerr,
      PetscScalar &gmaxHerr, PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
      PetscScalar &centerHerr) {
  // compute errors in thickness, eta=thickness^{(2n+2)/n}, volume, area
  
  PetscErrorCode  ierr;
  PetscScalar     **H, **HexactL;
  PetscScalar     Hexact, vol, area, domeH, volexact, areaexact, domeHexact;
  PetscScalar     Herr, avHerr, etaerr;

  PetscScalar     dummy, z, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  if (testname == 'L') {
    ierr = DAVecGetArray(grid.da2, vHexactL, &HexactL); CHKERRQ(ierr);
  }

  vol = 0; area = 0; domeH = 0;
  volexact = 0; areaexact = 0; domeHexact = 0;
  Herr = 0; avHerr=0; etaerr = 0;

  // area of grid square in square km:
  const PetscScalar   a = grid.p->dx * grid.p->dy * 1e-3 * 1e-3;
  const PetscScalar   m = (2*tgaIce.exponent()+2)/tgaIce.exponent();
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        area += a;
        vol += a * H[i][j] * 1e-3;
      }
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'A':
          exactA(r,&Hexact,&dummy);
          break;
        case 'B':
          exactB(grid.p->year*secpera,r,&Hexact,&dummy);
          break;
        case 'C':
          exactC(grid.p->year*secpera,r,&Hexact,&dummy);
          break;
        case 'D':
          exactD(grid.p->year*secpera,r,&Hexact,&dummy);
          break;
        case 'E':
          exactE(xx,yy,&Hexact,&dummy,&dummy1,&dummy2,&dummy3);
          break;
        case 'F':
          if (r > LforFG - 1.0) {  // outside of sheet
            Hexact=0.0;
          } else {
            r=PetscMax(r,1.0);
            z=0.0;
            bothexact(0.0,r,&z,1,0.0,
                      &Hexact,&dummy,&dummy5,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        case 'G':
          if (r > LforFG -1.0) {  // outside of sheet
            Hexact=0.0;
          } else {
            r=PetscMax(r,1.0);
            z=0.0;
            bothexact(grid.p->year*secpera,r,&z,1,ApforG,
                      &Hexact,&dummy,&dummy5,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        case 'H':
          exactH(f,grid.p->year*secpera,r,&Hexact,&dummy);
          break;
        case 'K':
          Hexact = 3000.0;
          break;
        case 'L':
          Hexact = HexactL[i][j];
          break;
        default:  SETERRQ(1,"test must be A, B, C, D, E, F, G, H, K, or L");
      }

      if (Hexact > 0) {
        areaexact += a;
        volexact += a * Hexact * 1e-3;
      }
      if (i == (grid.p->Mx - 1)/2 && j == (grid.p->My - 1)/2) {
        domeH = H[i][j];
        domeHexact = Hexact;
      }
      // compute maximum errors
      Herr = PetscMax(Herr,PetscAbsReal(H[i][j] - Hexact));
      etaerr = PetscMax(etaerr,PetscAbsReal(pow(H[i][j],m) - pow(Hexact,m)));
      // add to sums for average errors
      avHerr += PetscAbsReal(H[i][j] - Hexact);
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  if (testname == 'L') {
    ierr = DAVecRestoreArray(grid.da2, vHexactL, &HexactL); CHKERRQ(ierr);
  }
  
  // globalize (find errors over all processors) 
  PetscScalar gvol, garea, gdomeH;
  ierr = PetscGlobalSum(&volexact, &gvolexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeHexact, &gdomeHexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&areaexact, &gareaexact, grid.com); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&vol, &gvol, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  volerr = PetscAbsReal(gvol - gvolexact);
  areaerr = PetscAbsReal(garea - gareaexact);

  ierr = PetscGlobalMax(&Herr, &gmaxHerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avHerr, &gavHerr, grid.com); CHKERRQ(ierr);
  gavHerr = gavHerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalMax(&etaerr, &gmaxetaerr, grid.com); CHKERRQ(ierr);
  
  ierr = PetscGlobalMax(&domeH, &gdomeH, grid.com); CHKERRQ(ierr);
  centerHerr = PetscAbsReal(gdomeH - gdomeHexact);
  
  return 0;
}


PetscErrorCode IceCompModel::computeBasalVelocityErrors(
      PetscScalar &exactmaxspeed,
      PetscScalar &gmaxvecerr, PetscScalar &gavvecerr,
      PetscScalar &gmaxuberr, PetscScalar &gmaxvberr) {

  PetscErrorCode ierr;
  PetscScalar    **H, **ub, **vb;
  PetscScalar    maxvecerr, avvecerr, maxuberr, maxvberr;
  PetscScalar    ubexact,vbexact, dummy1,dummy2,dummy3;
  
  if (testname != 'E')
    SETERRQ(1,"basal velocity errors only computable for test E\n");
    
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  maxvecerr = 0.0; avvecerr = 0.0; maxuberr = 0.0; maxvberr = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      if (H[i][j] > 0.0) {
        PetscScalar r,xx,yy;
        mapcoords(i,j,xx,yy,r);
        exactE(xx,yy,&dummy1,&dummy2,&dummy3,&ubexact,&vbexact); 
        // compute maximum errors
        const PetscScalar uberr = PetscAbsReal(ub[i][j] - ubexact);
        const PetscScalar vberr = PetscAbsReal(vb[i][j] - vbexact);
        maxuberr = PetscMax(maxuberr,uberr);
        maxvberr = PetscMax(maxvberr,vberr);
        const PetscScalar vecerr = sqrt(uberr*uberr + vberr*vberr);
        maxvecerr = PetscMax(maxvecerr,vecerr);
        avvecerr += vecerr;      
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxuberr, &gmaxuberr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxvberr, &gmaxvberr, grid.com); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.p->Mx*grid.p->My);
  
  const PetscScalar pi = 3.14159265358979;
  const PetscScalar xpeak = 450e3 * cos(25.0*(pi/180.0)),
                    ypeak = 450e3 * sin(25.0*(pi/180.0));
  exactE(xpeak,ypeak,&dummy1,&dummy2,&dummy3,&ubexact,&vbexact);
  exactmaxspeed = sqrt(ubexact*ubexact + vbexact*vbexact);
  return 0;
}


PetscErrorCode IceCompModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

  PetscErrorCode ierr;
  if (printPrototype == PETSC_TRUE) {
    if ((testname == 'F') || (testname == 'G') || (testname == 'K')) {
      ierr = verbPrintf(2,grid.com,
               "P         YEAR:      ivol    iarea meltfABS    thick0     temp0\n");
      ierr = verbPrintf(2,grid.com,
               "U        years  10^6_km^3 10^6_km^2  (none)         m         K\n");
    } else {
      ierr = verbPrintf(2,grid.com,
               "P         YEAR:      ivol    iarea    thick0\n");
      ierr = verbPrintf(2,grid.com,
               "U        years  10^6_km^3 10^6_km^2        m\n");
    }
  } else {
    if ((testname == 'F') || (testname == 'G') || (testname == 'K')) {
      if (tempAndAge == PETSC_TRUE) {
        ierr = verbPrintf(2,grid.com, "S %12.5f: %9.5f %8.4f %8.4f %9.3f %9.4f\n",
                       year, volume_kmcube/1.0e6,area_kmsquare/1.0e6,meltfrac,H0,T0); CHKERRQ(ierr);
      } else {
        ierr = verbPrintf(2,grid.com, "S %12.5f: %9.5f %8.4f   <same> %9.3f    <same>\n",
                       year, volume_kmcube/1.0e6,area_kmsquare/1.0e6,H0); CHKERRQ(ierr);
      }
    } else {
        ierr = verbPrintf(2,grid.com, "S %12.5f: %9.5f %8.4f %9.3f\n",
           year, volume_kmcube/1.0e6, area_kmsquare/1.0e6, H0); CHKERRQ(ierr);
    }
  }
  return 0;
}


PetscErrorCode IceCompModel::reportErrors() {
  // geometry errors to report (for all tests): 
  //    -- max thickness error
  //    -- average (at each grid point on whole grid) thickness error
  //    -- dome thickness error
  //    -- max (thickness)^(2n+2)/n error
  //    -- volume error
  //    -- area error
  // and temperature errors (for tests F & G):
  //    -- max T error over 3D domain of ice
  //    -- av T error over 3D domain of ice
  // and basal temperature errors (for tests F & G):
  //    -- max basal temp error
  //    -- average (at each grid point on whole grid) basal temp error
  //    -- dome basal temp error
  // and strain-heating (Sigma) errors (for tests F & G):
  //    -- max Sigma error over 3D domain of ice (in 10^-3 K a^-1)
  //    -- av Sigma error over 3D domain of ice (in 10^-3 K a^-1)
  // and surface velocity errors (for tests F & G):
  //    -- max |<us,vs> - <usex,vsex>| error
  //    -- av |<us,vs> - <usex,vsex>| error
  //    -- max ws error
  //    -- av ws error
  // and basal sliding errors (for test E):
  //    -- max ub error
  //    -- max vb error
  //    -- max |<ub,vb> - <ubexact,vbexact>| error
  //    -- av |<ub,vb> - <ubexact,vbexact>| error

  PetscErrorCode  ierr;
  ierr = verbPrintf(1,grid.com, "NUMERICAL ERRORS evaluated at final time (relative to exact solution):\n"); CHKERRQ(ierr);

  // geometry (thickness, vol) errors if appropriate
  if (testname != 'K') {
    PetscScalar volexact, areaexact, domeHexact, volerr, areaerr, maxHerr, avHerr,
                maxetaerr, centerHerr;
    ierr = computeGeometryErrors(volexact,areaexact,domeHexact,
                                 volerr,areaerr,maxHerr,avHerr,maxetaerr,centerHerr);
            CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
            "geometry  :    prcntVOL        maxH         avH   relmaxETA      centerH\n");
            CHKERRQ(ierr);
    const PetscScalar   m = (2*tgaIce.exponent()+2)/tgaIce.exponent();
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f%13.6f\n",
                      100*volerr/volexact, maxHerr, avHerr,
                      maxetaerr/pow(domeHexact,m), centerHerr); CHKERRQ(ierr);
  }

  // temperature errors if appropriate
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxTerr, avTerr, basemaxTerr, baseavTerr, basecenterTerr;
    ierr = computeTemperatureErrors(maxTerr, avTerr); CHKERRQ(ierr);
    ierr = computeBasalTemperatureErrors(basemaxTerr, baseavTerr, basecenterTerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "temp      :        maxT         avT    basemaxT     baseavT  basecenterT\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f%13.6f\n", 
                  maxTerr, avTerr, basemaxTerr, baseavTerr, basecenterTerr); CHKERRQ(ierr);
  } else if (testname == 'K') {
    PetscScalar maxTerr, avTerr, maxTberr, avTberr;
    ierr = computeIceBedrockTemperatureErrors(maxTerr, avTerr, maxTberr, avTberr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "temp      :        maxT         avT       maxTb        avTb\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n", 
                  maxTerr, avTerr, maxTberr, avTberr); CHKERRQ(ierr);
  }

  // Sigma errors if appropriate
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxSigerr, avSigerr;
    ierr = computeSigmaErrors(maxSigerr, avSigerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "Sigma (3D):      maxSig       avSig\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f\n", 
                  maxSigerr*secpera*1.0e3, avSigerr*secpera*1.0e3); CHKERRQ(ierr);
  }

  // surface velocity errors if exact values are available
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxUerr, avUerr, maxWerr, avWerr;
    ierr = computeSurfaceVelocityErrors(maxUerr, avUerr, maxWerr, avWerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "surf vels :     maxUvec      avUvec        maxW         avW\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n", 
                  maxUerr*secpera, avUerr*secpera, maxWerr*secpera, avWerr*secpera); CHKERRQ(ierr);
  }

  // basal velocity errors if appropriate
  if (testname == 'E') {
    PetscScalar exactmaxspeed, maxvecerr, avvecerr, maxuberr, maxvberr;
    ierr = computeBasalVelocityErrors(exactmaxspeed,
                          maxvecerr,avvecerr,maxuberr,maxvberr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "base vels :  maxvector   avvector  prcntavvec     maxub     maxvb\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %11.4f%11.5f%12.5f%10.4f%10.4f\n", 
                  maxvecerr*secpera, avvecerr*secpera, 
                  (avvecerr/exactmaxspeed)*100.0,
                  maxuberr*secpera, maxvberr*secpera); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceCompModel::additionalAtStartTimestep() {
  PetscErrorCode    ierr;
  
  ierr = verbPrintf(5,grid.com,
               "additionalAtStartTimestep() in IceCompModel entered with test %c",
               testname); CHKERRQ(ierr);

  if (exactOnly == PETSC_TRUE)
    dt_force = maxdt;

  // these have no changing boundary conditions or comp sources:
  if (strchr("ABEKL",testname) != NULL) 
    return 0;

  switch (testname) {
    case 'C':
    case 'D':
    case 'H':
      ierr = getCompSourcesTestCDH();
      break;
    case 'F':
    case 'G':
      ierr = getCompSourcesTestFG();  // see iCMthermo.cc
      break;
    default:
      SETERRQ(1,"only tests CDHFG have comp source update at start time step\n");
      break;
  }

  return 0;
}


PetscErrorCode IceCompModel::additionalAtEndTimestep() {
  PetscErrorCode    ierr;
  
  ierr = verbPrintf(5,grid.com,
               "additionalAtEndTimestep() in IceCompModel entered with test %c",testname);
               CHKERRQ(ierr);

  // do nothing at the end of the time step unless the user has asked for the 
  // exact solution to overwrite the numerical solution
  if (exactOnly == PETSC_FALSE)  
    return 0;

  // because user want exact solution, fill gridded values from exact formulas
  switch (testname) {
    case 'A':
    case 'B':
    case 'C':
    case 'D':
    case 'H':
      ierr = fillSolnTestABCDH();
      break;
    case 'E':
      ierr = fillSolnTestE();
      break;
    case 'F':
    case 'G':
      ierr = fillSolnTestFG();  // see iCMthermo.cc
      break;
    case 'K':
      ierr = fillSolnTestK();  // see iCMthermo.cc
      break;
    case 'L':
      ierr = fillSolnTestL();
      break;
    default:
      SETERRQ(1,"unknown testname in IceCompModel");
      break;
  }

  return 0;
}

