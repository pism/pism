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
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../exact/exactTestIJ.h"
#include "iceExactSSAModel.hh"

const PetscScalar IceExactSSAModel::m_schoof = 10; // (pure number)
const PetscScalar IceExactSSAModel::L_schoof = 40e3; // meters
const PetscScalar IceExactSSAModel::aspect_schoof = 0.05; // (pure)
const PetscScalar IceExactSSAModel::H0_schoof = aspect_schoof * L_schoof; // = 2000 m THICKNESS
const PetscScalar IceExactSSAModel::B_schoof = 3.7e8; // Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big?
const PetscScalar IceExactSSAModel::p_schoof = 4.0/3.0; // = 1 + 1/n


IceExactSSAModel::IceExactSSAModel(IceGrid &g, IceType &i, char mytest)
  : IceModel(g,i) {
  test = mytest;
}


PetscErrorCode IceExactSSAModel::initFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      inFileSet, bifFileSet;
  char            inFile[PETSC_MAX_PATH_LEN];

  // does user want to turn off actual numerical evolution and simply report the
  //    exact solution ?
  ierr = PetscOptionsHasName(PETSC_NULL, "-eo", &exactOnly); CHKERRQ(ierr);

  // input file not allowed
  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bifFileSet); CHKERRQ(ierr);
  if ((inFileSet == PETSC_TRUE) || (bifFileSet == PETSC_TRUE)) {
    SETERRQ(2,"PISM input file not allowed for initialization of IceExactSSAModel");
  }
  
  ierr = verbPrintf(2,grid.com,"initializing Test I ... \n"); CHKERRQ(ierr);

  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  initialized_p = PETSC_TRUE;
  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);
  
  // set up 1000km by 240km by 3000m grid; note 240km = 6L where L = L_schoof
  // try -Mx 101 -My 25 for 10km x 10km grid
  ierr = grid.rescale(500e3, 120e3, 3000); CHKERRQ(ierr);

  // fill in temperature and age; not critical I think
  const PetscScalar T0 = 263.15;  // completely arbitrary
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);
  ierr = VecSet(vtau, 0.0); CHKERRQ(ierr);  // age, not yield stress

  // so make sure we are using plastic till
  if (createBasal_done == PETSC_TRUE) delete basal;
  basal = new PlasticBasalType;
  createBasal_done = PETSC_TRUE;
  
  ierr = taucSet(); CHKERRQ(ierr);  // now fill vtauc with values for Schoof manufactured solution
  
  return 0;
}


PetscErrorCode IceExactSSAModel::taucSet() {
  PetscErrorCode ierr;
  PetscScalar **tauc;

  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      const PetscScalar jfrom0
        = static_cast<PetscScalar>(j) - static_cast<PetscScalar>(grid.p->My - 1)/2.0;
      const PetscScalar y = grid.p->dy * jfrom0;
      const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const PetscScalar f = ice.rho * grav * H0_schoof * tan(theta);
      tauc[i][j] = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateAndBoundaryVels() {
  PetscErrorCode ierr;
  const PetscScalar    Mx = grid.p->Mx, My = grid.p->My;
  PetscScalar    **uvbar[2], **mask, **h, **bed;
  
  // set initial velocities in shelf for iteration
  ierr = VecSet(vubar,0.0 / secpera); CHKERRQ(ierr);
  ierr = VecSet(vvbar,0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);
  ierr = VecSet(vMask,MASK_DRAGGING); CHKERRQ(ierr);
  ierr = VecSet(vH,H0_schoof); CHKERRQ(ierr);

  // set h, bed everywhere
  // on edges y = +- 3 L_schoof, set velocity and make mask=SHEET
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
      bool edge = ( (j == 0) || (j == My-1) );
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

  // Communicate so that we can differentiate surface and set boundary conditions
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactSSAModel::reportErrors() {
  PetscErrorCode  ierr;
  const PetscScalar    Mx = grid.p->Mx, My = grid.p->My;
  PetscScalar exactmaxu, maxvecerr = 0.0, avvecerr = 0.0, 
              avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
              gmaxuerr = 0.0, gmaxverr = 0.0;
  PetscScalar **u, **v;

  ierr = verbPrintf(1,grid.com, 
          "NUMERICAL ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2, uexact, vexact;
      const PetscScalar 
              ifrom0 = static_cast<PetscScalar>(i) - static_cast<PetscScalar>(Mx - 1)/2.0,
              jfrom0 = static_cast<PetscScalar>(j) - static_cast<PetscScalar>(My - 1)/2.0;
      const PetscScalar myx = grid.p->dx * ifrom0, myy = grid.p->dy * jfrom0;
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &junk1, &junk2, &uexact, &vexact); 
      // compute maximum errors
      const PetscScalar uerr = PetscAbsReal(u[i][j] - uexact);
      const PetscScalar verr = PetscAbsReal(v[i][j] - vexact);
      avuerr = avuerr + uerr;      
      avverr = avverr + verr;      
      maxuerr = PetscMax(maxuerr,uerr);
      maxverr = PetscMax(maxverr,verr);
      const PetscScalar vecerr = sqrt(uerr * uerr + verr * verr);
      maxvecerr = PetscMax(maxvecerr,vecerr);
      avvecerr = avvecerr + vecerr;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
 
  // also get max of u
  PetscScalar junk1, junk2, junk3;
  exactI(m_schoof, 0.0, 0.0, &junk1, &junk2, &exactmaxu, &junk3);
  
  ierr = PetscGlobalMax(&maxuerr, &gmaxuerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxverr, &gmaxverr, grid.com); CHKERRQ(ierr);
  
  ierr = PetscGlobalSum(&avuerr, &gavuerr, grid.com); CHKERRQ(ierr);
  gavuerr = gavuerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
  gavverr = gavverr/(grid.p->Mx*grid.p->My);

  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.p->Mx*grid.p->My);

  ierr = verbPrintf(1,grid.com, 
     "      maxvector   avvector  prcntavvec      maxu      maxv       avu       avv\n");
     CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, 
            "    %11.4f%11.5f%12.5f%10.4f%10.4f%10.4f%10.4f\n", 
            gmaxvecerr*secpera, gavvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
            gmaxuerr*secpera, gmaxverr*secpera, gavuerr*secpera, gavverr*secpera); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, 
     "(exact maximum of u is %11.4f (m/a))\n",exactmaxu*secpera); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::diagnosticRun() {
  PetscErrorCode  ierr;
  PetscInt        pause_time;
  PetscTruth      pause_p;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, &pause_p); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running Test I ...\n"); CHKERRQ(ierr);
  if (exactOnly == PETSC_TRUE) {
    ierr=verbPrintf(2,grid.com,"  EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION\n"); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2,grid.com,
  "$$$$$       YEAR (+     STEP[$]):     VOL    AREA    MELTF     THICK0     TEMP0\n");
  CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep!

  ierr = setInitStateAndBoundaryVels(); CHKERRQ(ierr);

  // set flags, parameters affecting solve of stream equations
  useSSAVelocity = PETSC_TRUE;
  computeSurfGradInwardSSA = PETSC_TRUE;  // so periodic grid works although h(-Lx,y) != h(Lx,y)
  useConstantNuForSSA = PETSC_FALSE;
  useConstantHardnessForSSA = PETSC_TRUE;
  constantHardnessForSSA = B_schoof;
  ssaMaxIterations = 500;  
  setSSAEpsilon(0.0);  // don't use this lower bound

  if (exactOnly == PETSC_TRUE) { // just fill with exact solution
    PetscScalar **u, **v, **bed;
    ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        PetscScalar junk;
        const PetscScalar 
            ifrom0 = static_cast<PetscScalar>(i)
                       - static_cast<PetscScalar>(grid.p->Mx - 1)/2.0,
            jfrom0 = static_cast<PetscScalar>(j)
                       - static_cast<PetscScalar>(grid.p->My - 1)/2.0;
        const PetscScalar myx = grid.p->dx * ifrom0, myy = grid.p->dy * jfrom0;
        exactI(m_schoof, myx, myy, &bed[i][j], &junk, &u[i][j], &v[i][j]); 
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(4,grid.com,
                      "  [using Schoof regularization constant = %10.5e (m^2/s^2)\n",
                      PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof));
    CHKERRQ(ierr);
    ierr = basal->printInfo(4, grid.com); CHKERRQ(ierr);

    // solve model equations
    ierr = velocitySSA(); CHKERRQ(ierr);
  }

  // report on result of computation (i.e. to standard out and to viewers)
  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  ierr = summary(true,true); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);
  
  if (pause_p == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }
  return 0;
}
