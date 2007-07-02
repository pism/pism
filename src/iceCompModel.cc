// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cmath>
#include <cstring>
#include <petscda.h>

#include "exact/exactTestsABCDE.h"
#include "exact/exactTestsFG.h" 
#include "exact/exactTestH.h" 

#include "iceCompModel.hh"

// boundary conditions for tests F, G (same as EISMINT II Experiment F)
PetscScalar IceCompModel::Ggeo = 0.042;
PetscScalar IceCompModel::ST = 1.67e-5;
PetscScalar IceCompModel::Tmin = 223.15;  // K
PetscScalar IceCompModel::LforFG = 750000; // m
PetscScalar IceCompModel::ApforG = 200; // m
PetscScalar IceCompModel::ablationRateOutside = 0.02; // m/a

IceCompModel::IceCompModel(IceGrid &g, ThermoGlenArrIce &i)
  : IceModel(g, i), tgaIce(i) {
  
  // Override some defaults from parent class
  setEnhancementFactor(1.0);
  setThermalBedrock(PETSC_FALSE);
  setUseMacayealVelocity(PETSC_FALSE);
  setIsDrySimulation(PETSC_TRUE);
  setIncludeBMRinContinuity(PETSC_FALSE);

  f = tgaIce.rho / bedrock.rho;
  
  // now make bedrock have same material properties as ice
  // (note Mbz=1 also, by default, but want ice/rock interface to be pure ice)
  bedrock.rho = tgaIce.rho;
  bedrock.c_p = tgaIce.c_p;
  bedrock.k = tgaIce.k;

  // defaults for verification
  setTest('A');
  setExactOnly(PETSC_FALSE);
  compVecsCreated = PETSC_FALSE;
  compViewersCreated = PETSC_FALSE;
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
}


PetscErrorCode IceCompModel::destroyCompViewers() {
  PetscErrorCode ierr;
  if (SigmaCompView != PETSC_NULL) { ierr = PetscViewerDestroy(SigmaCompView); CHKERRQ(ierr); }
  if (compSigmaMapView != PETSC_NULL) { ierr = PetscViewerDestroy(compSigmaMapView); CHKERRQ(ierr); }
  return 0;
}


PetscErrorCode IceCompModel::createCompVecs() {
  PetscErrorCode ierr = DACreateLocalVector(grid.da3, &vSigmaComp); CHKERRQ(ierr);
  compVecsCreated = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceCompModel::destroyCompVecs() {
  PetscErrorCode  ierr = VecDestroy(vSigmaComp); CHKERRQ(ierr);
  return 0;
}


void IceCompModel::setTest(char c) {
  testname = c;
}


PetscErrorCode IceCompModel::setExactOnly(PetscTruth eo) {
  exactOnly = eo;
  return 0;
}


PetscErrorCode IceCompModel::setFromOptions() {
  PetscErrorCode ierr;
  char temptestname[20];
  char temp;

  ierr = IceModel::setFromOptions();  CHKERRQ(ierr);

  /* This option determines the single character name of a verification test
  ("-test B", for example). */
  ierr = PetscOptionsGetString(PETSC_NULL, "-test", temptestname, 1, &testchosen); CHKERRQ(ierr);
  if (testchosen == PETSC_TRUE) {
    temp = temptestname[0];
    if ((temp >= 'a') && (temp <= 'z'))
      temp += 'A'-'a';  // capitalize if lower
    if ((temp < 'A') || (temp > 'H'))
      SETERRQ(1,"IceCompModel ERROR: desired test NOT IMPLEMENTED\n");
  } else {
    temp = 'A';
  }
  setTest(temp);
  
  /* This switch turns off actual numerical evolution and simply reports the
     exact solution. */
  ierr = PetscOptionsHasName(PETSC_NULL, "-eo", &exactOnly); CHKERRQ(ierr);
  return 0;
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
  if ((testname == 'A') || (testname == 'E'))
    setOceanKill(PETSC_TRUE);

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, "continuing from input file %s; using Test %c conditions during run ...\n",inFile,testname);  CHKERRQ(ierr);
    grid.p->year=startYear; // some exact solutions have "absolute time"
  } else {
    ierr = verbPrintf(2,grid.com, "initializing Test %c ...\n",testname);  CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
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
        // use 1800km by 1800km by 4000m rectangular domain
        ierr = grid.rescale(900e3, 900e3, 4000); CHKERRQ(ierr);
        break;
      case 'H':
        // use 1500km by 1500km by 4000m rectangular domain
        ierr = grid.rescale(1500e3, 1500e3, 4000); CHKERRQ(ierr);
        break;
      default:  SETERRQ(1,"IceCompModel ERROR : desired test not implemented\n");
    }
    ierr = createVecs(); CHKERRQ(ierr);

    // none use Goldsby-Kohlstedt or do age calc
    ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS);
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    // all have no uplift at start
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);
    ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);

    grid.p->year=startYear; // some exact solutions have "absolute time"
    switch (testname) {
      case 'A':
      case 'B':
      case 'C':
      case 'D':
      case 'E':
      case 'H':
        ierr = initTestISO(); CHKERRQ(ierr);
        break;
      case 'F':
      case 'G':
        ierr = initTestFG(); CHKERRQ(ierr);
        break;
      default:  SETERRQ(1,"verify ERROR : desired test not implemented\n");
    }
  }

  ierr = afterInitHook(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::afterInitHook() {
  PetscErrorCode ierr;

  ierr = IceModel::afterInitHook(); CHKERRQ(ierr);
  
  ierr = createCompVecs(); CHKERRQ(ierr);
  ierr = createCompViewers();

  return 0;
}


PetscErrorCode IceCompModel::createCompViewers() {
  PetscErrorCode ierr;

  // must be called after IceModel::createViewers because diagnostic needs to be filled
  if (((testname=='F') || (testname=='G')) && (strchr(diagnostic, 'P') != NULL)) {
    ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, "Sigma_C (comPensatory heat) at kd",
             PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &SigmaCompView); CHKERRQ(ierr);
  } else SigmaCompView = PETSC_NULL;
  
  // take over SigmaMapView to show only strain heating and not sum Sigma + Sigma_C
  if (SigmaMapView != PETSC_NULL) {
    ierr = PetscViewerDestroy(SigmaMapView); CHKERRQ(ierr);
    SigmaMapView = PETSC_NULL;
    ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, "Sigma (strain heat) at kd",
             PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &compSigmaMapView); CHKERRQ(ierr);
  } else compSigmaMapView = PETSC_NULL;
   
  compViewersCreated = PETSC_TRUE;
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


// reimplement IceModel::basal
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


PetscErrorCode IceCompModel::initTestISO() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, **H, **accum, **mask, dummy1, dummy2, dummy3;
  const PetscScalar LforAE = 750e3; // m

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = isothermalFlux_A_softness;
  T0 = -tgaIce.Q() / (gasConst_R * log(A0/tgaIce.A()));
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);
  
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  setMuSliding(0.0);  // note reimplementation of basal()

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


PetscErrorCode IceCompModel::updateTestISO() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **accum, dummy, dummy1, dummy2, dummy3;

  // before flow step, set accumulation from exact values;
  // if exactOnly then also compute H here
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (exactOnly==PETSC_TRUE) {  // update H and accumulation
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
          case 'E':
            exactE(xx,yy,&H[i][j],&accum[i][j],&dummy,&dummy2,&dummy3);
            break;
          case 'H':
            exactH(f,grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
            break;
          default:  SETERRQ(1,"test must be A, B, C, D, E, or H");
        }
      } else { // real verification; only update accumulation
        switch (testname) {
          case 'A':
            accum[i][j] = 0.3/secpera; // no need to waste time ...
            break;
          case 'B':
            accum[i][j] = 0.0; // ditto ...
            break;
          case 'C':
            exactC(grid.p->year*secpera,r,&dummy,&accum[i][j]);
            break;
          case 'D':
            exactD(grid.p->year*secpera,r,&dummy,&accum[i][j]);
            break;
          case 'E':
            exactE(xx,yy,&dummy,&accum[i][j],&dummy1,&dummy2,&dummy3);
            break;
          case 'H':
            exactH(f,grid.p->year*secpera,r,&dummy,&accum[i][j]);
            break;
          default:  SETERRQ(1,"test must be A, B, C, D, E, or H");
        }
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);

  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
    if (testname == 'H') {
      ierr = VecCopy(vH,vh); CHKERRQ(ierr);
      ierr = VecScale(vh,1-f); CHKERRQ(ierr);
      ierr = VecCopy(vH,vbed); CHKERRQ(ierr);
      ierr = VecScale(vbed,-f); CHKERRQ(ierr);
    } else {
      ierr = VecCopy(vH,vh); CHKERRQ(ierr);
      ierr = VecSet(vbed,0.0); CHKERRQ(ierr);
    }
  }

  return 0;
}


PetscErrorCode IceCompModel::initTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.p->Mz;
  PetscScalar     **H, **accum, **Ts, ***T, ***Tb;
  PetscScalar     *z, *dummy1, *dummy2, *dummy3, *dummy4;
  z=new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  ierr = VecSet(vbed, 0); CHKERRQ(ierr);
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      Ts[i][j] = Tmin + ST * r;
      if (r > LforFG - 1.0) { // if (essentially) outside of sheet
        H[i][j]=0.0;
        accum[i][j] = -ablationRateOutside/secpera;
        for (PetscInt k=0; k<Mz; k++)
          T[i][j][k]=Ts[i][j];
      } else {
        r=PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
           bothexact(0.0,r,z,Mz,0.0,
                     &H[i][j],&accum[i][j],T[i][j],dummy1,dummy2,dummy3,dummy4);
        } else {
           bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                     &H[i][j],&accum[i][j],T[i][j],dummy1,dummy2,dummy3,dummy4);
        }
      }
      // fill with basal temp increased by geothermal flux
      for (PetscInt k=0; k<grid.p->Mbz; k++)
        Tb[i][j][k] = T[i][j][0] + bedrock.k * (grid.p->Mbz - k - 1) * grid.p->dz * Ggeo;
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);

  ierr = VecCopy(vH,vh); CHKERRQ(ierr);

  delete [] z;  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  return 0;
}


PetscErrorCode IceCompModel::updateTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.p->Mz;
  PetscScalar     **H, **accum, ***T, ***w, ***u, ***v, ***Sigma, ***SigmaC;
  PetscScalar     Ts, dummy0;
  PetscScalar     *z, *dummy1, *dummy2, *dummy3, *dummy4, *Uradial;
  z=new PetscScalar[Mz];  Uradial=new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  // before temperature and flow step, set Sigma_c and accumulation from
  // exact values for test F; also add SigmaC to Sigma for temperature computation
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);

  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  }

  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (r > LforFG - 1.0) {  // outside of sheet
        accum[i][j] = -ablationRateOutside/secpera;
        for (PetscInt k=0; k<Mz; k++)  SigmaC[i][j][k]=0.0;
        if (exactOnly == PETSC_TRUE) {
          H[i][j] = 0.0;
          Ts = Tmin + ST * r;
          for (PetscInt k=0; k<Mz; k++) {
            T[i][j][k]=Ts;
            u[i][j][k]=0.0;
            v[i][j][k]=0.0;
            w[i][j][k]=0.0;
            Sigma[i][j][k]=0.0;
          }
        }
      } else {
        r=PetscMax(r,1.0); // avoid singularity at origin
        if (exactOnly==PETSC_TRUE) {
          if (testname == 'F') {
            bothexact(0.0,r,z,Mz,0.0,
                      &H[i][j],&accum[i][j],T[i][j],Uradial,w[i][j],
                      Sigma[i][j],SigmaC[i][j]);
          } else {
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &H[i][j],&accum[i][j],T[i][j],Uradial,w[i][j],
                      Sigma[i][j],SigmaC[i][j]);
          }
          for (PetscInt k=0; k<Mz; k++) {
            u[i][j][k]=Uradial[k]*(xx/r);
            v[i][j][k]=Uradial[k]*(yy/r);
          }
        } else { // actually do update with compensatory Sigma; not exactOnly
          if (testname == 'F') {
            bothexact(0.0,r,z,Mz,0.0,
                      &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC[i][j]);
          } else {
             bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC[i][j]);
          }
          for (PetscInt k=0; k<Mz; k++)
            Sigma[i][j][k] += SigmaC[i][j][k]; // for rest of model, Sigma is both terms
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vSigma, INSERT_VALUES, vSigma); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vSigma, INSERT_VALUES, vSigma); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vSigmaComp, INSERT_VALUES, vSigmaComp); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vSigmaComp, INSERT_VALUES, vSigmaComp); CHKERRQ(ierr);

  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);
    ierr = VecCopy(vH,vh); CHKERRQ(ierr);
  }

  delete [] z;  delete [] Uradial;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  return 0;
}


PetscErrorCode IceCompModel::updateCompViewers() {
  PetscErrorCode ierr;
  PetscScalar    ***SigmaC, **sigmaC, ***Sigma, **sigma;

  ierr = updateViewers();  CHKERRQ(ierr);
  
  if (SigmaCompView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &sigmaC); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        sigmaC[i][j] = SigmaC[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &sigmaC); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, SigmaCompView); CHKERRQ(ierr);
  }
  
  // now redraw Sigma view to be just the strain-heating, not the sum of strain-heating 
  //      and compensatory
  if (compSigmaMapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &sigma); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        sigma[i][j] = Sigma[i][j][kd] - SigmaC[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &sigma); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, compSigmaMapView); CHKERRQ(ierr);
  }
  
  return 0;
}


PetscErrorCode IceCompModel::computeGeometryErrors(
      PetscScalar &gvolexact, PetscScalar &gareaexact, PetscScalar &gdomeHexact,
      PetscScalar &volerr, PetscScalar &areaerr,
      PetscScalar &gmaxHerr, PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
      PetscScalar &centerHerr) {
  // compute errors in thickness, thickness^{(2n+2)/n}, volume, area
  
  PetscErrorCode  ierr;
  PetscScalar     **H;
  PetscScalar     Hexact, vol, area, domeH, volexact, areaexact, domeHexact;
  PetscScalar     Herr, avHerr, etaerr;

  PetscScalar     dummy, z, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
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
        default:  SETERRQ(1,"test must be A, B, C, D, E, F, G, or H");
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


PetscErrorCode IceCompModel::computeTemperatureErrors(PetscScalar &gmaxTerr, PetscScalar &gavTerr) {

  PetscErrorCode ierr;
  PetscScalar    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  PetscScalar    **H, ***T;
  const PetscInt Mz = grid.p->Mz;
  
  PetscScalar   *z, *dummy1, *dummy2, *dummy3, *dummy4, *Tex;
  PetscScalar   junk0, junk1;
  z = new PetscScalar[Mz];  Tex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];  dummy4 = new PetscScalar[Mz];
  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        switch (testname) {
          case 'F':
            bothexact(0.0,r,z,Mz,0.0,
                      &junk0,&junk1,Tex,dummy1,dummy2,dummy3,dummy4);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &junk0,&junk1,Tex,dummy1,dummy2,dummy3,dummy4);
            break;
          default:  SETERRQ(1,"temperature errors only computable for tests F and G\n");
        }
       const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; k++) {  // only eval error if below num surface
          const PetscScalar Terr = PetscAbs(T[i][j][k] - Tex[k]);
          maxTerr = PetscMax(maxTerr,Terr);
          avcount += 1.0;
          avTerr += Terr;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

  delete [] z;  delete [] Tex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  
  ierr = PetscGlobalMax(&maxTerr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PetscGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeBasalTemperatureErrors(
      PetscScalar &gmaxTerr, PetscScalar &gavTerr, PetscScalar &centerTerr) {

  PetscErrorCode  ierr;
  PetscScalar     ***T, domeT, domeTexact, Terr, avTerr;

  PetscScalar     dummy, z, Texact, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  domeT=0; domeTexact = 0; Terr=0; avTerr=0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'F':
          if (r > LforFG - 1.0) {  // outside of sheet
            Texact=Tmin + ST * r;  // = Ts
          } else {
            r=PetscMax(r,1.0);
            z=0.0;
            bothexact(0.0,r,&z,1,0.0,
                      &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        case 'G':
          if (r > LforFG -1.0) {  // outside of sheet
            Texact=Tmin + ST * r;  // = Ts
          } else {
            r=PetscMax(r,1.0);
            z=0.0;
            bothexact(grid.p->year*secpera,r,&z,1,ApforG,
                      &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        default:  SETERRQ(1,"temperature errors only computable for tests F and G\n");
      }

      if (i == (grid.p->Mx - 1)/2 && j == (grid.p->My - 1)/2) {
        domeT = T[i][j][0];
        domeTexact = Texact;
      }
      // compute maximum errors
      Terr = PetscMax(Terr,PetscAbsReal(T[i][j][0] - Texact));
      // add to sums for average errors
      avTerr += PetscAbs(T[i][j][0] - Texact);
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  
  PetscScalar gdomeT, gdomeTexact;

  ierr = PetscGlobalMax(&Terr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalMax(&domeT, &gdomeT, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeTexact, &gdomeTexact, grid.com); CHKERRQ(ierr);  
  centerTerr = PetscAbsReal(gdomeT - gdomeTexact);
  
  return 0;
}


PetscErrorCode IceCompModel::computeSigmaErrors(PetscScalar &gmaxSigmaerr, PetscScalar &gavSigmaerr) {

  PetscErrorCode ierr;
  PetscScalar    maxSigerr = 0.0, avSigerr = 0.0, avcount = 0.0;
  PetscScalar    **H, ***Sig;
  const PetscInt Mz = grid.p->Mz;
  
  PetscScalar   *z, *dummy1, *dummy2, *dummy3, *Sigex, *SigCompex;
  PetscScalar   junk0, junk1;
  z = new PetscScalar[Mz];  Sigex = new PetscScalar[Mz];  SigCompex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];
  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sig); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        switch (testname) {
          case 'F':
            bothexact(0.0,r,z,Mz,0.0,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,SigCompex);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,SigCompex);
            break;
          default:  SETERRQ(1,"strain-heating (Sigma) errors only computable for tests F and G\n");
        }
        // verbPrintf(1,grid.com,"%e|",Sigex[3]);
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; k++) {  // only eval error if below num surface
          const PetscScalar actualSignum  = Sig[i][j][k] - SigCompex[k];
          const PetscScalar Sigerr = PetscAbs(actualSignum - Sigex[k]);
          maxSigerr = PetscMax(maxSigerr,Sigerr);
          avcount += 1.0;
          avSigerr += Sigerr;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sig); CHKERRQ(ierr);

  delete [] z;  delete [] Sigex;  delete [] SigCompex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;
  
  ierr = PetscGlobalMax(&maxSigerr, &gmaxSigmaerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avSigerr, &gavSigmaerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PetscGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavSigmaerr = gavSigmaerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeSurfaceVelocityErrors(
        PetscScalar &gmaxUerr, PetscScalar &gavUerr,
        PetscScalar &gmaxWerr, PetscScalar &gavWerr) {

  PetscErrorCode ierr;
  PetscScalar    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;
  PetscScalar    **H, **unum, **vnum, **wnum;

  ierr = getSurfaceValuesOf3D(vu, vWork2d[0]); CHKERRQ(ierr); // = numerical surface val of u
  ierr = getSurfaceValuesOf3D(vv, vWork2d[1]); CHKERRQ(ierr); // = numerical surface val of v
  ierr = getSurfaceValuesOf3D(vw, vWork2d[2]); CHKERRQ(ierr); // = numerical surface val of w

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &unum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &vnum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &wnum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                               // and not at central singularity
        PetscScalar radialUex,wex;
        PetscScalar dummy0,dummy1,dummy2,dummy3,dummy4;
        PetscScalar z = H[i][j];
        switch (testname) {
          case 'F':
            bothexact(0.0,r,&z,1,0.0,
                      &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,&z,1,ApforG,
                      &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);
            break;
          default:  SETERRQ(1,"surface velocity errors only computed for tests F and G\n");
        }
        // verbPrintf(1,grid.com,"[%f|%f]",radialUex,wex);
        const PetscScalar uex = (xx/r) * radialUex;
        const PetscScalar vex = (yy/r) * radialUex;
        const PetscScalar Uerr = sqrt(PetscSqr(unum[i][j] - uex) + PetscSqr(vnum[i][j] - vex));
        maxUerr = PetscMax(maxUerr,Uerr);
        avUerr += Uerr;
        const PetscScalar Werr = PetscAbs(wnum[i][j] - wex);
        maxWerr = PetscMax(maxWerr,Werr);
        avWerr += Werr;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &unum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &vnum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &wnum); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxUerr, &gmaxUerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxWerr, &gmaxWerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avUerr, &gavUerr, grid.com); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalSum(&avWerr, &gavWerr, grid.com); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.p->Mx*grid.p->My);
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

  // geometry (thickness, vol) errors; area not reported as meaningless as computed
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

  // temperature errors if appropriate
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxTerr, avTerr, basemaxTerr, baseavTerr, basecenterTerr;
    ierr = computeTemperatureErrors(maxTerr, avTerr); CHKERRQ(ierr);
    ierr = computeBasalTemperatureErrors(basemaxTerr, baseavTerr, basecenterTerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "temp      :        maxT         avT    basemaxT     baseavT  basecenterT\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f%13.6f\n", 
                  maxTerr, avTerr, basemaxTerr, baseavTerr, basecenterTerr); CHKERRQ(ierr);
  }

  // Sigma errors if appropriate
  if ((testname == 'F') || (testname == 'G')) {
    PetscScalar maxSigerr, avSigerr;
    ierr = computeSigmaErrors(maxSigerr, avSigerr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
       "Sigma (3D):         max          av\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f\n", 
                  maxSigerr*secpera*1.0e3, avSigerr*secpera*1.0e3); CHKERRQ(ierr);
  }

  // surface velocity errors if appropriate
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


PetscErrorCode IceCompModel::dumpToFile_Matlab(const char *fname) {
  PetscErrorCode  ierr;
  PetscViewer  viewer;
  char b[PETSC_MAX_PATH_LEN];  // for basename
  char mcompf[PETSC_MAX_PATH_LEN];  // Matlab format for compensatory

  ierr = IceModel::dumpToFile_Matlab(fname);  CHKERRQ(ierr);
  
  if ((testname == 'F') || (testname == 'G')) {
    ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
    strcpy(mcompf, b);
    strcat(mcompf, "_comp.m");

    ierr = verbPrintf(2,grid.com, "\n(also dumping Sigma_C to `%s'; also corrects Sigma)", mcompf); CHKERRQ(ierr);

    ierr = PetscViewerASCIIOpen(grid.com, mcompf, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"\n\n\n\ndisp('iceCompModel output:')\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nclear SigmaCkd\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

    PetscScalar     ***SigmaC, **SigmaC2;
    ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &SigmaC2); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        SigmaC2[i][j] = SigmaC[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &SigmaC2); CHKERRQ(ierr);
    ierr=PetscObjectSetName((PetscObject) g2,"SigmaCkd"); CHKERRQ(ierr);
    ierr = LVecView(grid.da2, vWork2d[0],  g2, viewer); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nSigmaCkd = reshape(SigmaCkd,%d,%d);\n\n",
          grid.p->Mx,grid.p->My);  CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nSigmakd = Sigmakd - SigmaCkd; \n\n");  CHKERRQ(ierr);
 
    ierr = PetscViewerASCIIPrintf(viewer,"figure\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"mesh(x,y,Sigmakd), colormap cool\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"title('TRUE strain heating \\Sigma at z given by -kd')\n\n");
         CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"figure\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"mesh(x,y,SigmaCkd), colormap cool\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"title('compensatory heating \\Sigma_C at z given by -kd')\n\n");
         CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);
  }
  return 0;
}


PetscErrorCode IceCompModel::run() {
  PetscErrorCode  ierr;
  
  ierr = verbPrintf(2,grid.com, "running test %c ...\n", testname); CHKERRQ(ierr);
  if (exactOnly == PETSC_TRUE) {
    ierr=verbPrintf(2,grid.com,"  EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,grid.com,
      "$$$$       YEAR (+     STEP[N$]):     VOL    AREA MELTFabs     THICK0     TEMP0\n");
      CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"$$$$");  CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  tempskipCountDown = 0;
  ierr = summary(true,false); CHKERRQ(ierr);  // report starting state

  dtTempAge = 0.0;
  // main loop for time evolution
  for (PetscScalar year = startYear; year < endYear; year += dt/secpera) {
    dt_force = -1.0;
    maxdt_temporary = -1.0;
    ierr = additionalAtStartTimestep(); CHKERRQ(ierr);  // might set dt_force,maxdt_temp

    // compute bed deformation, which only depends on current thickness and bed elevation
    if (doBedDef == PETSC_TRUE) {
      ierr = bedDefStepIfNeeded(); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

    // always do vertically-average velocity calculation; only update velocities at depth if
    // needed for temp and age calculation
//    bool tempAgeStep = (    (exactOnly == PETSC_FALSE)
//                         && (doTemp == PETSC_TRUE)
//                         && (tempskipCountDown == 0) 
//                        && ((testname == 'F') || (testname =='G')) );
//    ierr = velocity(tempAgeStep); CHKERRQ(ierr);
    bool updateAtDepth = (    (exactOnly == PETSC_FALSE)
                          && (tempskipCountDown == 0));
    ierr = velocity(updateAtDepth); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, updateAtDepth ? "v" : "V" ); CHKERRQ(ierr);

    // adapt time step using velocities and diffusivity, ..., just computed
    if (exactOnly == PETSC_TRUE)
      dt_force = maxdt;
    bool useCFLforTempAgeEqntoGetTimestep = 
              ((doTemp == PETSC_TRUE) && ((testname == 'F') || (testname == 'G')));
    ierr = determineTimeStep(useCFLforTempAgeEqntoGetTimestep); CHKERRQ(ierr);
    dtTempAge += dt;
    grid.p->year += dt / secpera;  // adopt it
    // IceModel::dt,dtTempAge,grid.p->year are now set correctly according to
    //    mass-continuity-eqn-diffusivity criteria, CFL criteria, and other 
    //    criteria from derived class additionalAtStartTimestep(), and from 
    //    "-tempskip" mechanism

    // ierr = PetscPrintf(PETSC_COMM_SELF,
    //           "\n[rank=%d, it=%d, year=%f, dt=%f]", grid.rank, it, year, dt/secpera);
    //        CHKERRQ(ierr);
    

    switch (testname) {
      case 'A':
      case 'B':
      case 'C':
      case 'D':
      case 'E':
      case 'H':
        ierr = updateTestISO(); CHKERRQ(ierr);
        break;
      case 'F':
      case 'G':
        ierr = updateTestFG(); CHKERRQ(ierr);
        break;
      default: SETERRQ(1, "test must be A, B, C, D, E, F, G, or H\n");
    }
    
    bool tempAgeStep = ( (doTemp == PETSC_TRUE)
                         && ((testname == 'F') || (testname =='G'))
                         && updateAtDepth );
    if (tempAgeStep) {
      // note temps are allowed to go above pressure melting in verify
      allowAboveMelting = PETSC_TRUE;
      ierr = temperatureAgeStep(); CHKERRQ(ierr);
      dtTempAge = 0.0;
      ierr = verbPrintf(2,grid.com, "t"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
    if ((exactOnly == PETSC_FALSE) && (doMassConserve == PETSC_TRUE)) {
      ierr = massBalExplicitStep(); CHKERRQ(ierr);
      if ((doTempSkip == PETSC_TRUE) && (tempskipCountDown > 0))
        tempskipCountDown--;
      ierr = verbPrintf(2,grid.com, "f"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

    ierr = summary(tempAgeStep,false); CHKERRQ(ierr);
    ierr = updateCompViewers(); CHKERRQ(ierr);
    
    ierr = additionalAtEndTimestep(); CHKERRQ(ierr);
    if (endOfTimeStepHook() != 0) break;
  }  //  for loop

  return 0;
}
