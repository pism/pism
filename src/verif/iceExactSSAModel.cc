// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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
#include "exactTestsIJ.h"
#include "iceExactSSAModel.hh"

const PetscScalar IceExactSSAModel::m_schoof = 10; // (pure number)
const PetscScalar IceExactSSAModel::L_schoof = 40e3; // meters
const PetscScalar IceExactSSAModel::aspect_schoof = 0.05; // (pure)
const PetscScalar IceExactSSAModel::H0_schoof = aspect_schoof * L_schoof; 
                                       // = 2000 m THICKNESS
const PetscScalar IceExactSSAModel::B_schoof = 3.7e8; // Pa s^{1/3}; hardness 
                                       // given on p. 239 of Schoof; why so big?
const PetscScalar IceExactSSAModel::p_schoof = 4.0/3.0; // = 1 + 1/n

const PetscScalar IceExactSSAModel::LforJ = 300.0e3; // 300 km half-width


IceExactSSAModel::IceExactSSAModel(IceGrid &g, IceType *i, char mytest)
  : IceModel(g,i) {
  test = mytest;
}


PetscErrorCode IceExactSSAModel::initFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      inFileSet, bifFileSet;
  char            inFile[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(2,grid.com,"initializing Test %c ... \n",test); CHKERRQ(ierr);

  // input file not allowed
  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bifFileSet); CHKERRQ(ierr);
  if ((inFileSet == PETSC_TRUE) || (bifFileSet == PETSC_TRUE)) {
    SETERRQ(2,"PISM input file not allowed for initialization of IceExactSSAModel");
  }
  
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
    
  // fill in temperature and age; not critical
  const PetscScalar T0 = 263.15;  // completely arbitrary
  ierr = VecSet(vTs, T0); CHKERRQ(ierr);
  ierr = T3.setToConstant(T0); CHKERRQ(ierr);
  ierr = Tb3.setToConstant(T0); CHKERRQ(ierr);

  ierr = tau3.setToConstant(0.0); CHKERRQ(ierr);  // age (not yield stress)

  // set initial velocities in shelf (for start of iteration)
  ierr = VecSet(vubar,0.0); CHKERRQ(ierr);
  ierr = VecSet(vvbar,0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);
  // clear 3D and basal velocities too
  ierr = u3.setToConstant(0.0); CHKERRQ(ierr);
  ierr = v3.setToConstant(0.0); CHKERRQ(ierr);
  ierr = w3.setToConstant(0.0); CHKERRQ(ierr);
  ierr = VecSet(vub,0.0); CHKERRQ(ierr);
  ierr = VecSet(vvb,0.0); CHKERRQ(ierr);
  
  useSSAVelocity = PETSC_TRUE;
  ssaMaxIterations = 500;  
  useConstantNuForSSA = PETSC_FALSE;
  doPlasticTill = PETSC_TRUE;  // correct for I, irrelevant for J
  doSuperpose = PETSC_FALSE;
  
  // do rest of base class initialization before final initialization for I and J (below)
  initialized_p = PETSC_TRUE;
  ierr = IceModel::initFromOptions(PETSC_FALSE); CHKERRQ(ierr);
  // in preparation for rescale_and_set_zlevels() below:
  ierr = determineSpacingTypeFromOptions(); CHKERRQ(ierr);  

  switch (test) {
    case 'I': {
      // set up X km by 240 km by 3000 m grid; note 240km = 6L where L = L_schoof
      // find dimension X so that pixels are square but for some unknown reason 
      //    there are convergence problems when x dimension Lx is significantly
      //    smaller than y dimension Ly
      // note Mx,My may have been set by options, but not others
      const PetscScalar testI_Ly = 120.0e3,
        testI_Lx = PetscMax(60.0e3, ((grid.Mx - 1) / 2) * (2.0 * testI_Ly / (grid.My - 1)) );
      ierr = grid.rescale_and_set_zlevels(testI_Lx, testI_Ly, 3000.0); CHKERRQ(ierr);
      ierr = afterInitHook(); CHKERRQ(ierr);
      // fill vtauc with values for Schoof solution:
      ierr = taucSetI(); CHKERRQ(ierr);
      // set up remaining stuff:
      ierr = setInitStateAndBoundaryVelsI(); CHKERRQ(ierr);
      // set flags, parameters affecting solve of stream equations
      // so periodic grid works although h(-Lx,y) != h(Lx,y):
      computeSurfGradInwardSSA = PETSC_TRUE;
      useConstantHardnessForSSA = PETSC_TRUE;
      constantHardnessForSSA = B_schoof;
      ssaEpsilon = 0.0;  // don't use this lower bound
      
      ierr = verbPrintf(3,grid.com,
      "IceExactSSAModel::initFromOptions:regularizingVelocitySchoof = %10.5e\n",
      regularizingVelocitySchoof); CHKERRQ(ierr);
      }
      break;
    case 'J':
      // first allocate space for nu (which will be calculated from formula)
      ierr = VecDuplicateVecs(vh, 2, &vNuForJ); CHKERRQ(ierr);
      // we set a flag because the grid is truely periodic
      ierr = grid.rescale_and_set_zlevels(LforJ, LforJ, 1000, PETSC_TRUE, PETSC_TRUE);
         CHKERRQ(ierr);
      ierr = afterInitHook(); CHKERRQ(ierr);
      ierr = setInitStateJ(); CHKERRQ(ierr);
      isDrySimulation = PETSC_FALSE;
      leaveNuAloneSSA = PETSC_TRUE; // will use already-computed nu instead of updating
      computeSurfGradInwardSSA = PETSC_FALSE;
      useConstantHardnessForSSA = PETSC_FALSE;
      ssaEpsilon = 0.0;  // don't use this lower bound
      break;
    default:
      SETERRQ(1,"Only tests I and J currently supported in IceExactSSAModel.");
  }

  // Communicate so that we can differentiate surface, and to set boundary conditions
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceExactSSAModel::taucSetI() {
  PetscErrorCode ierr;
  PetscScalar **tauc;

  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt jfrom0 = j - (grid.My - 1)/2;
      const PetscScalar y = grid.dy * jfrom0;
      const PetscScalar theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
      const PetscScalar f = ice->rho * grav * H0_schoof * tan(theta);
      tauc[i][j] = f * pow(PetscAbs(y / L_schoof), m_schoof);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::setInitStateAndBoundaryVelsI() {
  PetscErrorCode ierr;
  PetscScalar    **uvbar[2], **mask, **h, **bed;
  
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
      const PetscInt ifrom0 = i - (grid.Mx - 1)/2,
                     jfrom0 = j - (grid.My - 1)/2;
      const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
      // eval exact solution; will only use exact vels if at edge
      exactI(m_schoof, myx, myy, &(bed[i][j]), &junk, &myu, &myv); 
      h[i][j] = bed[i][j] + H0_schoof;
      bool edge = ( (j == 0) || (j == grid.My - 1) );
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


PetscErrorCode IceExactSSAModel::setInitStateJ() {
  PetscErrorCode ierr;
  PetscScalar    **H, **h, **mask, **uvbar[2];

  ierr = VecSet(vbed,-5000.0); CHKERRQ(ierr); // assures shelf is floating
  ierr = VecSet(vMask,MASK_FLOATING); CHKERRQ(ierr);

  /* use Ritz et al (2001) value of 30 MPa yr for typical vertically-averaged viscosity */
  const PetscScalar nu0 = 30.0 * 1.0e6 * secpera; /* = 9.45e14 Pa s */
  const PetscScalar H0 = 500.0;       /* 500 m typical thickness */
  ierr = VecSet(vNuForJ[0], nu0 * H0); CHKERRQ(ierr);
  ierr = VecSet(vNuForJ[1], nu0 * H0); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar junk1, myu, myv;
      const PetscInt ifrom0 = i - (grid.Mx)/2,
                     jfrom0 = j - (grid.My)/2;
      const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
      // set H,h on regular grid
      ierr = exactJ(myx, myy, &H[i][j], &junk1, &myu, &myv); CHKERRQ(ierr);
      h[i][j] = (1.0 - ice->rho / ocean.rho) * H[i][j];
      // special case at center point: here we indirectly set ubar,vbar 
      // at (i,j) by marking this grid point as SHEET and setting staggered-grid
      // version of ubar approriately; the average done in assembleSSARhs()
      // for SHEET points puts our value of velocity onto regular grid point (i,j)
      // and makes ubar=myu, vbar=myv there
      if ( (i == (grid.Mx)/2) && (j == (grid.My)/2) ) {
        mask[i][j] = MASK_SHEET;
        uvbar[0][i-1][j] = myu;   uvbar[0][i][j] = myu;
        uvbar[1][i][j-1] = myv;   uvbar[1][i][j] = myv;
        
        ierr = u3.needAccessToVals(); CHKERRQ(ierr);
        ierr = v3.needAccessToVals(); CHKERRQ(ierr);
        ierr = u3.setToConstantColumn(i,j,myu); CHKERRQ(ierr);
        ierr = v3.setToConstantColumn(i,j,myv); CHKERRQ(ierr);
        ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
        ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
      }
    }
  }  
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::reportErrors() {
  PetscErrorCode  ierr;
  const PetscInt    Mx = grid.Mx, My = grid.My;
  PetscScalar exactmaxu, maxvecerr = 0.0, avvecerr = 0.0, 
              avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
              gmaxuerr = 0.0, gmaxverr = 0.0;
  PetscScalar **u, **v;

  if (doPseudoPlasticTill == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com, 
          "WARNING: numerical errors not valid for pseudo-plastic till\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com, 
          "NUMERICAL ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2, uexact, vexact;
      const PetscInt ifrom0 = (test == 'I') ? i - (Mx - 1)/2 : i - (Mx)/2,
                     jfrom0 = (test == 'I') ? j - (My - 1)/2 : j - (My)/2;
      const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
      // eval exact solution
      if (test == 'I') {
        exactI(m_schoof, myx, myy, &junk1, &junk2, &uexact, &vexact); 
      } else if (test == 'J') {
        exactJ(myx, myy, &junk1, &junk2, &uexact, &vexact);
      } else {
        SETERRQ(1,"only tests I and J have computable errors");
      }
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
   
  ierr = PetscGlobalMax(&maxuerr, &gmaxuerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxverr, &gmaxverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avuerr, &gavuerr, grid.com); CHKERRQ(ierr);
  gavuerr = gavuerr/(grid.Mx*grid.My);
  ierr = PetscGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
  gavverr = gavverr/(grid.Mx*grid.My);
  ierr = PetscGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);

  if (test == 'I') {
    PetscScalar junk1, junk2, junk3;
    exactI(m_schoof, 0.0, 0.0, &junk1, &junk2, &exactmaxu, &junk3);
    ierr = verbPrintf(1,grid.com, 
       "velocity  :  maxvector   prcntavvec      maxu       avu\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
            "           %11.4f%13.5f%10.4f%10.4f\n", 
            gmaxvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
            gmaxuerr*secpera, gavuerr*secpera); CHKERRQ(ierr);
  } else if (test == 'J') {
    // following from "pismv -test J -Mx 601 -My 601 -Mz 3 -verbose -eo"
    exactmaxu = 181.366 / secpera;  
    ierr = verbPrintf(1,grid.com, 
       "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
       CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
            "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n", 
            gmaxvecerr*secpera, (gavvecerr/exactmaxu)*100.0,
            gmaxuerr*secpera, gmaxverr*secpera, gavuerr*secpera, 
            gavverr*secpera); CHKERRQ(ierr);
  }


  if (test == 'I') {
    // also print max of u
    ierr = verbPrintf(3,grid.com, 
       "(exact maximum of u is %11.4f (m/a))\n",exactmaxu*secpera); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceExactSSAModel::fillFromExactSolution() {
  PetscErrorCode  ierr;
  PetscScalar **ubar, **vbar;

  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar junk1, junk2;
      const PetscScalar 
          ifrom0 = static_cast<PetscScalar>(i)
                     - static_cast<PetscScalar>(grid.Mx - 1)/2.0,
          jfrom0 = static_cast<PetscScalar>(j)
                     - static_cast<PetscScalar>(grid.My - 1)/2.0;
      const PetscScalar myx = grid.dx * ifrom0, myy = grid.dy * jfrom0;
      if (test == 'I') {
        exactI(m_schoof, myx, myy, &junk1, &junk2, &ubar[i][j], &vbar[i][j]); 
      } else if (test == 'J') {
        ierr = exactJ(myx, myy, &junk1, &junk2, &ubar[i][j], &vbar[i][j]); CHKERRQ(ierr);
      } else {
        SETERRQ(1,"only tests I and J supported in IceExactSSAModel");
      }
      for (PetscInt k=0; k<grid.Mz; k++) {
        ierr = u3.setToConstantColumn(i,j,ubar[i][j]); CHKERRQ(ierr);
        ierr = v3.setToConstantColumn(i,j,vbar[i][j]); CHKERRQ(ierr);
        // leave w alone for now; will be filled by call to broadcastSSAVelocity()
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceExactSSAModel::diagnosticRun() {
  PetscErrorCode  ierr;

  ierr = verbPrintf(2,grid.com, "running Test %c ...\n", test); CHKERRQ(ierr);

  // does user want to turn off actual numerical evolution (and save exact solution)?
  ierr = PetscOptionsHasName(PETSC_NULL, "-eo", &exactOnly); CHKERRQ(ierr);
  if (exactOnly == PETSC_TRUE) {
    ierr=verbPrintf(2,grid.com,"  EXACT SOLUTION ONLY, NO NUMERICAL SOLUTION\n"); CHKERRQ(ierr);
  }

  // print out some stats about input state
  ierr = summaryPrintLine(PETSC_TRUE,PETSC_TRUE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
           CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  skipCountDown = 0;

  if (exactOnly == PETSC_TRUE) {
    // just fill with exact solution, including exact 3D hor velocities
    ierr = fillFromExactSolution(); CHKERRQ(ierr);
  } else { // numerically solve ice shelf/stream equations
    PetscInt numiter;
    ierr = initSSA(); CHKERRQ(ierr);
    if (test == 'I') {
      ierr = velocitySSA(&numiter); CHKERRQ(ierr);
    } else if (test == 'J') {
      // use locally allocated space for (computed) nu:
      ierr = velocitySSA(vNuForJ,&numiter); CHKERRQ(ierr); 
      ierr = VecDestroyVecs(vNuForJ, 2); CHKERRQ(ierr); // immediately de-allocate
    }
    // fill in 3D velocities (u,v,w)
    ierr = broadcastSSAVelocity(true); CHKERRQ(ierr);
  }

  // finally update w
  ierr = u3.beginGhostComm(); CHKERRQ(ierr);
  ierr = u3.endGhostComm(); CHKERRQ(ierr);
  ierr = v3.beginGhostComm(); CHKERRQ(ierr);
  ierr = v3.endGhostComm(); CHKERRQ(ierr);
  ierr = vertVelocityFromIncompressibility(); CHKERRQ(ierr);

  // report on result of computation (i.e. to standard out and to viewers)
  ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  ierr = summary(true,true); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);
  PetscInt    pause_time = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, PETSC_NULL); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }
  return 0;
}

