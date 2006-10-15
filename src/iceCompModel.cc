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

#include <cmath>
#include <cstring>
#include <petscda.h>

#include "exactTestsBCD.h"
#include "exactTestsFG.h" 
// WHAT IS GOING ON?  WHY WON'T IT LINK IF I INCLUDE THE HEADER FILES?

#include "iceCompModel.hh"

// boundary conditions for tests F, G (i.e. EISMINT II Experiment F)
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
  setMuSliding(0.0);
  setThermalBedrock(PETSC_FALSE);
  setUseMacayealVelocity(PETSC_FALSE);
  setIsDrySimulation(PETSC_TRUE);
  setIsothermalFlux(PETSC_FALSE);

  f = tgaIce.rho / bedrock.rho;

  // Defaults specific to this model
  setTest('F');
  setExactOnly(PETSC_FALSE);
}


IceCompModel::~IceCompModel() {
  destroyCompViewers();
  destroyCompVecs();
}


PetscErrorCode IceCompModel::createCompVecs() {
  PetscErrorCode ierr = DACreateLocalVector(grid.da3, &vSigmaComp); CHKERRQ(ierr);
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
  
  /* This option determines the single character name of a verification test:
  "-ts B", for example. */
  ierr = PetscOptionsGetString(PETSC_NULL, "-test", temptestname, 1, &testchosen); CHKERRQ(ierr);
  if (testchosen == PETSC_TRUE) {
    temp = temptestname[0];
    if ((temp >= 'a') && (temp <= 'z'))   temp += 'A'-'a';  // capitalize if lower
    switch (temp) {
 // list implemented tests here:
      case 'B':
      case 'C':
      case 'D':
      case 'F':
      case 'G':
        setTest(temp);
        break;
      case 'H':
        if ((doBedDef == PETSC_FALSE) || (doBedIso == PETSC_FALSE)) {
          ierr = PetscPrintf(grid.com, 
                "[verify WARNING: Test H should be run with option  -bed_def_iso  for the reported errors to be correct.]\n");
                CHKERRQ(ierr);
        }
        setTest(temp);
        break;
      default:  SETERRQ(1,"verify ERROR: desired test NOT IMPLEMENTED\n");
    }
  }
  
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

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, "continuing; using Test %c conditions ...\n",testname);  CHKERRQ(ierr);
    grid.p->year=startYear; // some exact solutions have "absolute time"
  } else {
    ierr = PetscPrintf(grid.com, "initializing Test %c..\n",testname);  CHKERRQ(ierr);
    ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);

    // none use Goldsby-Kohlstedt or do age calc
    ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS);
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    // all have no uplift at start
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);

    grid.p->year=startYear; // some exact solutions have "absolute time"
    switch (testname) {
      case 'B':
        // use 2400km by 2400km by 4000m rectangular domain
        ierr = grid.rescale(1200e3, 1200e3, 4000); CHKERRQ(ierr);
        ierr = initTestBCDH(); CHKERRQ(ierr);
        break;
      case 'C':
      case 'D':
        // use 2000km by 2000km by 4000m rectangular domain
        ierr = grid.rescale(1000e3, 1000e3, 4000); CHKERRQ(ierr);
        ierr = initTestBCDH(); CHKERRQ(ierr);
        break;
      case 'F':
      case 'G':
        // use 1800km by 1800km by 4000m rectangular domain
        ierr = grid.rescale(900e3, 900e3, 4000); CHKERRQ(ierr);
        ierr = initTestFG(); CHKERRQ(ierr);
        break;
      case 'H':
        // use 1500km by 1500km by 4000m rectangular domain
        ierr = grid.rescale(1500e3, 1500e3, 4000); CHKERRQ(ierr);
        ierr = initTestBCDH(); CHKERRQ(ierr);
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
   
  return 0;
}


PetscErrorCode IceCompModel::destroyCompViewers() {
  PetscErrorCode ierr;
  if (SigmaCompView != PETSC_NULL) { ierr = PetscViewerDestroy(SigmaCompView); CHKERRQ(ierr); }
  if (compSigmaMapView != PETSC_NULL) { ierr = PetscViewerDestroy(compSigmaMapView); CHKERRQ(ierr); }
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


int IceCompModel::exactH(double t, double r, double &H, double &M) {
  // just a wrapper for exactC_iso to form test H

  // f should be already set
  const double n = 3.0;
  // t0 = (beta/Gamma) * pow((2n+1)/((n+1)(1-f)),n) * (pow(R0,n+1)/pow(H0,2n+1)) when beta=2;
  double t0 = (15208.0 / pow(1-f,n)) * secpera;  // 40033 years  (for test C with isostasy)
  double lambda, alpha, beta, Rmargin;
  const double H0 = 3600.0, R0=750000.0;

  if (t < t0) { // t <= t0: version of test C
    lambda=5.0;
    alpha=-1.0;  // alpha=(2-(n+1)*lambda)/(5*n+3)
    beta=2.0;  // beta=(1+(2*n+1)*lambda)/(5*n+3)
  } else { // t >= t0: version of test B
    const double t0post = (t0/2.0) * (1.0/18.0);  // reset t and t0 
    t = t - t0 + t0post; // reset to Halfar w. f
    t0 = t0post;
    lambda=0.0;
    alpha=1.0/9.0;  // alpha=(2-(n+1)*lambda)/(5*n+3)=1/9
    beta=1.0/18.0;  // beta=(1+(2*n+1)*lambda)/(5*n+3)=1/18
  }

  Rmargin = R0 * pow(t/t0,beta);
  if (r < Rmargin)
    H = H0 * pow(t/t0,-alpha) * pow(  1.0-pow( pow(t/t0,-beta)*(r/R0), (n+1)/n ),  n/(2*n+1)  );
  else
    H = 0.0;

  if (t > 0.1*secpera)
    M = (lambda/t) * H;
  else {  // when less than 0.1 year, avoid division by time
    Rmargin = R0 * pow(0.1*secpera/t0,beta);
    if (r < Rmargin)
      M = lambda * H0 / t0;  // constant value in disc of Rmargin radius
    else
      M = 0.0; 
  }
  
  return 0;
}


PetscErrorCode IceCompModel::initTestBCDH() {
  PetscErrorCode  ierr;
  PetscScalar     A0, T0, **H, **accum;

  // compute T so that A0 = A(T) = Acold exp(-Qcold/(R T))  (i.e. for ThermoGlenArrIce);
  // set all temps to this constant
  A0 = isothermalFlux_A_softness;
  T0 = -tgaIce.Q() / (tgaIce.gasConst_R * log(A0/tgaIce.A()));
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);

  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (testname == 'B')
        exactB(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
      else if (testname == 'C')
        exactC(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
      else if (testname == 'D')
        exactD(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
      else if (testname == 'H') {
        exactH(grid.p->year*secpera,r,H[i][j],accum[i][j]);
      }
      else SETERRQ(1,"test must be B, C, D, or H");
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
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


PetscErrorCode IceCompModel::updateTestBCDH() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **accum, dummy;

  // before flow step, set accumulation from exact values; if exactOnly then also compute H here
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (exactOnly==PETSC_TRUE) {  // update H and accumulation
        if (testname == 'B')
          exactB(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
        else if (testname == 'C')
          exactC(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
        else if (testname == 'D')
          exactD(grid.p->year*secpera,r,&H[i][j],&accum[i][j]);
        else if (testname == 'H')
          exactH(grid.p->year*secpera,r,H[i][j],accum[i][j]);
        else SETERRQ(1,"test must be B, C, D, or H");
      } else {
        if (testname == 'B')
          accum[i][j] = 0.0;  // no need to waste time calling exactB() for this
        else if (testname == 'C')
          exactC(grid.p->year*secpera,r,&dummy,&accum[i][j]);
        else if (testname == 'D')
          exactD(grid.p->year*secpera,r,&dummy,&accum[i][j]);
        else if (testname == 'H')
          exactH(grid.p->year*secpera,r,dummy,accum[i][j]);
        else SETERRQ(1,"test must be B, C, D, or H");
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
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
  PetscScalar     **H, **accum, **Ts, ***T;
  PetscScalar     *z, *dummy1, *dummy2, *dummy3, *dummy4;
  z=new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  globalMinTemp=Tmin;

  ierr = VecSet(vbed, 0); CHKERRQ(ierr);
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);

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
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
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
  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
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


PetscErrorCode IceCompModel::reporterror() {

  PetscErrorCode  ierr;
  PetscScalar     **H, ***T;
  PetscScalar     Hexact, vol, area, domeH, domeT, volexact, areaexact, domeHexact, domeTexact;
  PetscScalar     Herr, avHerr, etaerr, Terr, avTerr;

  PetscScalar     dummy, *z, *Texact, *dummy1, *dummy2, *dummy3, *dummy4;
  const PetscInt  Mz=1;
  z=new PetscScalar[Mz];  Texact = new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  // geometry errors to report: 
  //    -- max thickness error
  //    -- average (at each grid point on whole grid) thickness error
  //    -- dome thickness error
  //    -- max (thickness)^(2n+2)/n error
  //    -- volume error
  //    -- area error
  // and temperature errors:
  //    -- max basal temp error [tests F&G]
  //    -- average (at each grid point on whole grid) basal temp error [tests F&G]
  //    -- dome basal temp error [tests F&G]

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  vol = 0; area = 0; domeH = 0; domeT=0; 
  volexact = 0; areaexact = 0; domeHexact = 0; domeTexact = 0;
  Herr = 0; avHerr=0; etaerr = 0; Terr=0; avTerr=0;

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
        case 'B':
          exactB(grid.p->year*secpera,r,&Hexact,&dummy);
          Texact[0] = T[i][j][0];  // no temp error to report
          break;
        case 'C':
          exactC(grid.p->year*secpera,r,&Hexact,&dummy);
          Texact[0] = T[i][j][0];  // no temp error to report
          break;
        case 'D':
          exactD(grid.p->year*secpera,r,&Hexact,&dummy);
          Texact[0] = T[i][j][0];  // no temp error to report
          break;
        case 'F':
          if (r > LforFG - 1.0) {  // outside of sheet
            Hexact=0.0;
            Texact[0]=Tmin + ST * r;  // = Ts
          } else {
            r=PetscMax(r,1.0);
            z[0]=0.0;
            bothexact(0.0,r,z,Mz,0.0,
                      &Hexact,&dummy,Texact,dummy1,dummy2,dummy3,dummy4);
          }
          break;
        case 'G':
          if (r > LforFG -1.0) {  // outside of sheet
            Hexact=0.0;
            Texact[0]=Tmin + ST * r;  // = Ts
          } else {
            r=PetscMax(r,1.0);
            z[0]=0.0;
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &Hexact,&dummy,Texact,dummy1,dummy2,dummy3,dummy4);
          }
          break;
        case 'H':
          exactH(grid.p->year*secpera,r,Hexact,dummy);
          Texact[0] = T[i][j][0];  // no temp error to report
          break;
      }

      if (Hexact > 0) {
        areaexact += a;
        volexact += a * Hexact * 1e-3;
      }
      if (i == (grid.p->Mx - 1)/2 && j == (grid.p->My - 1)/2) {
        domeH = H[i][j];
        domeHexact = Hexact;
        domeT = T[i][j][0];
        domeTexact = Texact[0];
      }
      // compute maximum errors
      Herr = PetscMax(Herr,PetscAbsReal(H[i][j] - Hexact));
      etaerr = PetscMax(etaerr,PetscAbsReal(pow(H[i][j],m) - pow(Hexact,m)));
      Terr = PetscMax(Terr,PetscAbsReal(T[i][j][0] - Texact[0]));
      // add to sums for average errors
      avHerr += PetscAbsReal(H[i][j] - Hexact);
      avTerr += PetscAbsReal(T[i][j][0] - Texact[0]);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  
  PetscScalar gvol, garea, gvolexact, gareaexact, 
              gdomeH, gdomeHexact, gdomeT, gdomeTexact,
              gHerr, gTerr, getaerr, gavHerr, gavTerr;

  ierr = PetscGlobalSum(&vol, &gvol, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volexact, &gvolexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&areaexact, &gareaexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&Herr, &gHerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&etaerr, &getaerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&Terr, &gTerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeH, &gdomeH, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeHexact, &gdomeHexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeT, &gdomeT, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeTexact, &gdomeTexact, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avHerr, &gavHerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  
  PetscScalar   volerr = PetscAbsReal(gvol - gvolexact);
  PetscScalar   areaerr = PetscAbsReal(garea - gareaexact);
  PetscScalar   domeHerr = PetscAbsReal(gdomeH - gdomeHexact);
  PetscScalar   domeTerr = PetscAbsReal(gdomeT - gdomeTexact);
  
  gavHerr = gavHerr/(grid.p->Mx*grid.p->My);
  gavTerr = gavTerr/(grid.p->Mx*grid.p->My);

  ierr = PetscPrintf(grid.com, "Actual ERRORS (evaluated at final time):\n"); CHKERRQ(ierr);
  // print geometry errors
  ierr = PetscPrintf(grid.com, 
     "geometry  :  prcntVOL  prcntAREA      maxH         avH   relmaxETA    domeH\n");
     CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com, "           %10.4f%11.4f%10.4f%12.6f%12.6f%9.4f\n",
                100*volerr/gvolexact, 100*areaerr/gareaexact, gHerr, gavHerr,
                getaerr/pow(gdomeHexact,m), domeHerr); CHKERRQ(ierr);
  // print temp errors
  ierr = PetscPrintf(grid.com, 
     "base temps:        maxT         avT      domeT\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com, "           %12.6f%12.6f%11.6f\n", 
                gTerr, gavTerr, domeTerr); CHKERRQ(ierr);

  delete [] z;  delete [] Texact;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;

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

    ierr = PetscPrintf(grid.com, "\n(also dumping Sigma_C to `%s'; also corrects Sigma)", mcompf); CHKERRQ(ierr);

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
  // 7/19/06: only tests B, C, D, F, G, H
  PetscErrorCode  ierr;
  
  ierr = PetscPrintf(grid.com, "running test %c ...\n", testname); CHKERRQ(ierr);
  if (exactOnly == PETSC_TRUE) {
    ierr=PetscPrintf(grid.com,"  EXACT SOLUTION ONLY, NO NUMERICS\n"); CHKERRQ(ierr);
  }
  ierr = initSounding(); CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com,
      "$$$$      YEAR (+   STEP):     VOL    AREA MELTFabs     THICK0     TEMP0\n");
      CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com,"$$$$");  CHKERRQ(ierr);
  ierr = summary(true,false); CHKERRQ(ierr);  // report starting state

  PetscScalar dt_temp = 0.0;
  PetscInt    it = 0;
  bool        tempAgeStep;
  
  // main loop
  for (PetscScalar year = startYear; year < endYear; year += dt/secpera, it++) {
    // compute bed deformation, which only depends on current thickness and bed elevation
    if (doBedDef == PETSC_TRUE) {
      ierr = bedDefStepIfNeeded(); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
    }

    // always do vertically-average velocity calculation; only update velocities at depth if
    // needed for temp and age calculation
    tempAgeStep = ( (exactOnly == PETSC_FALSE) && (doTemp == PETSC_TRUE) && 
                          (it % tempskip == 0) && ((testname == 'F') || (testname =='G')) );
    // if (doVelocity == PETSC_TRUE) // flag ignored
    ierr = velocity(tempAgeStep); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, tempAgeStep ? "v" : "V" ); CHKERRQ(ierr);

    // adapt time step using velocities just computed
    dt = PetscMin(maxdt, (endYear-year) * secpera);  // don't go past end; "propose" this
    if (doAdaptTimeStep == PETSC_TRUE) {
      if ( (exactOnly == PETSC_FALSE) && (doMassBal == PETSC_TRUE)) {
        ierr = adaptTimeStepDiffusivity();  CHKERRQ(ierr);
      }
      if ( (exactOnly == PETSC_FALSE) && (doTemp == PETSC_TRUE) && 
               ((testname == 'F') || (testname =='G')) ) {
        ierr = adaptTimeStepCFL();  CHKERRQ(ierr);  // if tempskip > 1 then here dt is reduced
                                                    // by a factor of tempskip
      }
    }    
    // IceModel::dt is now set correctly according to mass-balance and CFL criteria  
    dt_temp += dt;
    grid.p->year += dt / secpera;  // adopt it

    switch (testname) {
      case 'B':
      case 'C':
      case 'D':
        ierr = updateTestBCDH(); CHKERRQ(ierr);
        break;
      case 'F':
      case 'G':
        ierr = updateTestFG(); CHKERRQ(ierr);
        break;
      case 'H':
        ierr = updateTestBCDH(); CHKERRQ(ierr);
        break;
      default:
        SETERRQ(1, "Test must be B, C, D, F, G, H\n");
    }
    
    if (tempAgeStep) {
      // note temps are allowed to go above pressure melting in verify
      ierr = temperatureStep(PETSC_TRUE, dt_temp); CHKERRQ(ierr);  // also does age
      dt_temp = 0.0;
      ierr = PetscPrintf(grid.com, "t"); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
    }
    
    if ((exactOnly == PETSC_FALSE) && (doMassBal == PETSC_TRUE)) {
      ierr = massBalExplicitStep(); CHKERRQ(ierr);
      ierr = PetscPrintf(grid.com, "f"); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
    }

    ierr = summary(tempAgeStep,false); CHKERRQ(ierr);

    ierr = updateCompViewers(); CHKERRQ(ierr);
  }  //  for loop

  return 0;
}
