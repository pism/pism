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

#include <cstring>
#include <cmath>
#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

#include "iceExactStreamModel.hh"
#include "exact/exactTestI.h"

const PetscScalar IceExactStreamModel::m_schoof = 10; // (pure number)
const PetscScalar IceExactStreamModel::L_schoof = 40e3; // meters
const PetscScalar IceExactStreamModel::aspect_schoof = 0.05; // (pure)
const PetscScalar IceExactStreamModel::H0_schoof = aspect_schoof * L_schoof; // = 2000 m THICKNESS
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
  PetscTruth      sometestchosen, inFileSet;
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
  
  ierr = verbPrintf(2,grid.com,"initializing Test I ... \n"); CHKERRQ(ierr);

  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);

  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  
  // set up 1000km by 240km by 3000m grid; note 240km = 6L where L = L_schoof
  // try -Mx 101 -My 25 for 10km x 10km grid
  ierr = grid.rescale(500e3, 120e3, 3000); CHKERRQ(ierr);
  
  ierr = fillinTemps();  CHKERRQ(ierr);

  ierr = afterInitHook(); CHKERRQ(ierr);  // note this sets basal to ViscousBasalType

  // make sure we are using plastic till
  if (createBasal_done == PETSC_TRUE) delete basal;
  basal = new PlasticBasalType;
  createBasal_done = PETSC_TRUE;
  ierr = taucSet(); CHKERRQ(ierr);  // now fill vtauc with values for Schoof manufactured solution
  
  return 0;
}


PetscErrorCode IceExactStreamModel::fillinTemps() {
  PetscErrorCode      ierr;
  const PetscScalar T0 = 263.15;  // completely arbitrary
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = VecSet(vT, T0); CHKERRQ(ierr);
  ierr = VecSet(vTb, T0); CHKERRQ(ierr);
  ierr = VecSet(vtau, 0.0); CHKERRQ(ierr);  // age, not yield stress
  return 0;
}


PetscErrorCode IceExactStreamModel::taucSet() {
  // compare IceCompModel::mapcoords()
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


// reimplement basal drag as plastic law [THIS MAY OR MAY NOT WORK!!]
#if 0
PetscScalar IceExactStreamModel::basalDragx(PetscScalar **beta, PetscScalar **tauc,
                                            PetscScalar **u, PetscScalar **v,
                                            PetscInt i, PetscInt j) const {
  //return taucGet(i,j) / sqrt(PetscSqr(plastic_regularize) + PetscSqr(u[i][j]) + PetscSqr(v[i][j]));
  return basal->drag(beta[i][j], tauc[i][j], u[i][j], v[i][j]);
}


// ditto
PetscScalar IceExactStreamModel::basalDragy(PetscScalar **beta, PetscScalar **tauc,
                                            PetscScalar **u, PetscScalar **v,
                                            PetscInt i, PetscInt j) const {
  // return taucGet(i,j) / sqrt(PetscSqr(plastic_regularize) + PetscSqr(u[i][j]) + PetscSqr(v[i][j]));
  return basal->drag(beta[i][j], tauc[i][j], u[i][j], v[i][j]);
}
#endif

PetscErrorCode IceExactStreamModel::setInitStateAndBoundaryVels() {
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
//      bool edge = (   (i == 0) || (i == 1) || (i == Mx-2) || (i == Mx-1)
//                   || (j == 0) || (j == 1) || (j == My-2) || (j == My-1) );
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

  // Communicate so that we can successfully differentiate surface and set boundary conditions
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactStreamModel::reportErrors() {
  PetscErrorCode  ierr;
  const PetscScalar    Mx = grid.p->Mx, My = grid.p->My;
  PetscScalar exactmaxu, maxvecerr = 0.0, avvecerr = 0.0, 
              avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
              gmaxuerr = 0.0, gmaxverr = 0.0;
  PetscScalar **u, **v;

  ierr = verbPrintf(2,grid.com, 
          "Actual ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);

//  ierr = verbPrintf(2,grid.com, "  xs = %d, xs+xm = %d, ys = %d, ys+ym = %d, Mx * My = %d\n",
//            grid.xs,grid.xs+grid.xm,grid.ys,grid.ys+grid.ym,grid.p->Mx*grid.p->My); CHKERRQ(ierr);
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
//  gavuerr = avuerr;
  gavuerr = gavuerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
//  gavuerr = avverr;
  gavverr = gavverr/(grid.p->Mx*grid.p->My);

  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
//  gavvecerr = avvecerr;
  gavvecerr = gavvecerr/(grid.p->Mx*grid.p->My);

  ierr = verbPrintf(2,grid.com, 
     "      maxvector   avvector  prcntavvec      maxu      maxv       avu       avv\n");
     CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, 
            "    %11.4f%11.5f%12.5f%10.4f%10.4f%10.4f%10.4f\n", 
            gmaxvecerr*secpera, gavvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
            gmaxuerr*secpera, gmaxverr*secpera, gavuerr*secpera, gavverr*secpera); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, 
     "(exact maximum of u is %11.4f (m/a))\n",exactmaxu*secpera); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactStreamModel::run() {
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

  ierr = fillinTemps(); CHKERRQ(ierr);
  ierr = setInitStateAndBoundaryVels(); CHKERRQ(ierr);

  // set flags, parameters affecting solve of stream equations
  useMacayealVelocity = PETSC_TRUE;
  computeSurfGradInwardMacAyeal = PETSC_TRUE;  // so periodic grid works even though
                                               // h(-Lx,y) != h(Lx,y)
  useConstantNuForMacAyeal = PETSC_FALSE;
  useConstantHardnessForMacAyeal = PETSC_TRUE;
  constantHardnessForMacAyeal = B_schoof;
  macayealMaxIterations = 500;  
  setMacayealEpsilon(0.0);  // don't use this lower bound
  // regularizingVelocitySchoof = 1.0 / secpera;  // 1 m/a is small velocity for ice stream?
  // regularizingLengthSchoof = 1000.0e3;         // (VELOCITY/LENGTH)^2  is very close to 10^-27

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

    ierr = velocityMacayeal(); CHKERRQ(ierr);
  }

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
