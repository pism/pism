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

#include "iceHEINOModel.hh"

IceHEINOModel::IceHEINOModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {
}


void IceHEINOModel::setExperName(char name) {
  expername = name;
}


char IceHEINOModel::getExperName() {
  return expername;
}


void IceHEINOModel::setflowlawNumber(PetscInt law) {
  flowlawNumber = law;
}


PetscInt IceHEINOModel::getflowlawNumber() {
  return flowlawNumber;
}


PetscErrorCode IceHEINOModel::setExperNameFromOptions() {
  PetscErrorCode      ierr;
  char                ismipexpername[20], runName[20];
  PetscTruth          ISMIPchosen, datprefixchosen, ismiprunchosen;

  /* This option chooses ISMIP; "-ismip H" is ISMIP-HEINO */
  ierr = PetscOptionsGetString(PETSC_NULL, "-ismip", ismipexpername, 1, &ISMIPchosen); CHKERRQ(ierr);

/*  this block might make sense later, but for now just set experiment name to HEINO:
  if (ISMIPchosen == PETSC_TRUE) {
    char temp = ismipexpername[0];
    if ((temp == 'H') || (temp == 'h')) // ISMIP-HEINO
        setExperName('0');
    else if ((temp == 'O') || (temp == 'o')) {
        setExperName('1');
        SETERRQ(1, "ISMIP-HOM NOT IMPLEMENTED!");
    }
    else if ((temp == 'P') || (temp == 'p')) {
        setExperName('2');
        SETERRQ(1, "ISMIP-POLICE NOT IMPLEMENTED!");
    }
  } else { // set a default: ISMIP HEINO
    setExperName('0');
  }
*/
  setExperName('0');
  
  /* if this option is set then no .dat files are produced for ISMIP-HEINO runs */  
  ierr = PetscOptionsHasName(PETSC_NULL, "-no_deliver", &ismipNoDeliver); CHKERRQ(ierr);
  /* if this option is set then 0.25 year time steps are not forced (but a deliverable
     is written at each integer year) */  
  ierr = PetscOptionsHasName(PETSC_NULL, "-allow_adapt", &ismipAllowAdapt); CHKERRQ(ierr);
  /* see ISMIP-HEINO documentation; options are [corresponding val of ismipHeinoRun]:
       -run ST  [0; default]
       -run T1  [1]
       -run T2  [2]
       -run B1  [3]
       -run B2  [4]
       -run S1  [5]
       -run S2  [6]
       -run S3  [7]    */
  ierr = PetscOptionsGetString(PETSC_NULL, "-datprefix", heinodatprefix, 20,
                               &datprefixchosen);  CHKERRQ(ierr);
  if (datprefixchosen != PETSC_TRUE) {
    strcpy(heinodatprefix,"PISM");
  }
  ierr = PetscOptionsGetString(PETSC_NULL, "-run", runName, 20, &ismiprunchosen);
            CHKERRQ(ierr);
  for (int i = 0; i<(int) strlen(runName); i++) {
    if ((runName[i] >= 'a') && (runName[i] <= 'z'))
      runName[i] += 'A'-'a';  // capitalize if lower
  }
  if (ismiprunchosen != PETSC_TRUE) {
    ismipHeinoRun = 0;
    strcpy(ismipRunName,"ST");
  } else {
    if      (strstr(runName,"ST") != NULL) {
      ismipHeinoRun = 0;
      strcpy(ismipRunName,"ST");
    } else if (strstr(runName,"T1") != NULL) {
      ismipHeinoRun = 1;
      strcpy(ismipRunName,"T1");
    } else if (strstr(runName,"T2") != NULL) {
      ismipHeinoRun = 2;
      strcpy(ismipRunName,"T2");
    } else if (strstr(runName,"B1") != NULL) {
      ismipHeinoRun = 3;
      strcpy(ismipRunName,"B1");
    } else if (strstr(runName,"B2") != NULL) {
      ismipHeinoRun = 4;
      strcpy(ismipRunName,"B2");
    } else if (strstr(runName,"S1") != NULL) {
      ismipHeinoRun = 5;
      strcpy(ismipRunName,"S1");
    } else if (strstr(runName,"S2") != NULL) {
      ismipHeinoRun = 6;
      strcpy(ismipRunName,"S2");
    } else if (strstr(runName,"S3") != NULL) {
      ismipHeinoRun = 7;        
      strcpy(ismipRunName,"S3");
    } else {
      SETERRQ(1,"ERROR setExperNameFromOptions(): known ISMIP-HEINO name not found");
    }
  }
  return 0;
}


PetscErrorCode IceHEINOModel::initFromOptions() {
  PetscErrorCode      ierr;
  char                inFile[PETSC_MAX_PATH_LEN];
  const PetscScalar   G_geothermal   = 0.042;             // J/m^2 s; geo. heat flux

  ierr = setExperNameFromOptions(); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, 
              "initializing ISMIP-HEINO, run %s ... \n", 
              ismipRunName); CHKERRQ(ierr);
    ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
//    grid.p->Mbz = 0; // overrides options
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);
    // if no inFile then starts with zero ice
    ierr = VecSet(vh, 0);
    ierr = VecSet(vH, 0);
    ierr = VecSet(vbed, 0);
    ierr = VecSet(vHmelt, 0.0);
    ierr = VecSet(vGhf, G_geothermal);
    ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS);
    // none use Goldsby-Kohlstedt
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    // all have no uplift at start
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);

    // next block checks if experiment is implemented
    // note height of grid must be great enough to handle max thickness
    switch (getExperName()) {
      case '0': // ISMIP-HEINO
        ierr = grid.rescale(2000e3, 2000e3, 7000); CHKERRQ(ierr);
        break;
      default:  
        SETERRQ(1,"ERROR: desired ISMIP experiment NOT IMPLEMENTED\n");
    }
  }

  if (ismipNoDeliver != PETSC_TRUE) {
    ierr = heinoCreateDat(); CHKERRQ(ierr);
  }

  ierr = applyDefaultsForExperiment(); CHKERRQ(ierr);

  ierr = initAccumTs(); CHKERRQ(ierr);
  if (inFileSet == PETSC_FALSE) {
    ierr = fillintemps(); CHKERRQ(ierr);
  }

  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = PetscPrintf(grid.com, "running ISMIP-HEINO, run %s ...\n",ismipRunName);
             CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceHEINOModel::applyDefaultsForExperiment() {
  PetscErrorCode   ierr;
  PetscScalar      **mask;

  setUseMacayealVelocity(PETSC_FALSE);
  setIsDrySimulation(PETSC_TRUE);
  setDoGrainSize(PETSC_FALSE);
  setEnhancementFactor(3.0);
//  setIncludeBMRinContinuity(PETSC_FALSE);
  setIncludeBMRinContinuity(PETSC_TRUE);
  setOceanKill(PETSC_TRUE);

  // make bedrock material properties into ice properties
  // (note Mbz=0 always, but want ice/rock interface segment to be all ice)
  bedrock.rho = ice.rho;
  bedrock.k = ice.k;
  bedrock.c_p = ice.c_p;
  
  if (getExperName() == '0') { // ISMIP-HEINO

    C_S = 500 / secpera;
    if (ismipHeinoRun == 5) C_S = 100 / secpera; // S1
    if (ismipHeinoRun == 6) C_S = 200 / secpera; // S2
    if (ismipHeinoRun == 7) C_S = 1000 / secpera; // S3
        
    ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar myx = -grid.p->Lx + grid.p->dx * i, 
                          myy = -grid.p->Ly + grid.p->dy * j;
        const PetscScalar r = sqrt(PetscSqr(myx) + PetscSqr(myy));
        if (r > (2000e3 - 1))
          mask[i][j] = MASK_FLOATING_OCEAN0;
        else
          mask[i][j] = MASK_SHEET;
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  } else {
    SETERRQ(1,"NO OTHER ISMIP EXPERIMENTS IMPLEMENTED");
  }

  return 0;
}


PetscErrorCode IceHEINOModel::initAccumTs() {
  PetscErrorCode    ierr;
  PetscScalar       **accum, **Ts;

  const PetscScalar  S_T = 2.5e-9 / 1e9; // 2.5e-9 K km^-3
  PetscScalar  T_min = 233.15;
  if (ismipHeinoRun == 1) T_min = 223.15; // T1
  if (ismipHeinoRun == 2) T_min = 243.15; // T2
    
  PetscScalar  b_min = 0.15 / secpera;
  PetscScalar  b_max = 0.3 / secpera;
  if (ismipHeinoRun == 3) { // B1
    b_min = 0.075 / secpera;
    b_max = 0.15 / secpera;
  }
  if (ismipHeinoRun == 4) { // B2
    b_min = 0.3 / secpera;
    b_max = 0.6 / secpera;
  }

  // now fill in accum and surface temp
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // d is distance from center of grid
      const PetscScalar d = sqrt( PetscSqr(-grid.p->Lx + grid.p->dx*i)
                                  + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
      // set accumulation
      accum[i][j] = b_min + ((b_max - b_min) / 2000e3) * d;
      // set surface temperature
      Ts[i][j] = T_min + S_T * d * d * d;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceHEINOModel::fillintemps() {
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



bool IceHEINOModel::inSoftSediment(const PetscScalar x, const PetscScalar y) {
  const bool inbigrect = (   (x >= 2300e3-1) && (x <= 3300e3+1)
                          && (y >= 1500e3-1) && (y <= 2500e3+1) ),
             insmallrect = (   (x >= 3300e3-1) && (x <= 4000e3+1)
                            && (y >= 1900e3-1) && (y <= 2100e3+1) );
  return (inbigrect || insmallrect);
}


// reimplement IceModel::basal(): location-dependent pressure-melting-temperature-activated
  // linear *or nonlinear* sliding law.  Returns positive coefficient C in the law
  //                U_b = <u_b,v_b> = - C grad h 
  // note: ignors mu
PetscScalar IceHEINOModel::basal(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {
  
  if (getExperName() == '0') { // ISMIP-HEINO
    //PetscErrorCode  ierr = PetscPrintf(grid.com, 
    //        "   [IceHEINOModel::basal called with:   x=%f, y=%f, H=%f, T=%f, alpha=%f]\n",
    //        x,y,H,T,alpha);  CHKERRQ(ierr);
    if (T + ice.beta_CC_grad * H > DEFAULT_MIN_TEMP_FOR_SLIDING) {
      // set coords according to ISMIP-HEINO convention
      const PetscScalar  xIH = x + grid.p->Lx,
                         yIH = y + grid.p->Ly;
      if (inSoftSediment(xIH,yIH)) {
        // see applyDefaultsForExperiment() for setting C_S
        return C_S * H;
      } else {
        // "hard rock"; C_R = 10^5 a^-1; alpha = |grad h|
        return (1e5 /secpera) * H * alpha * alpha;
      }
    } else
      return 0.0;  // not at pressure melting
  } else {
    SETERRQ(1,"experiment name invalid in IceHEINOModel::basal()");
    return 0.0;
  }
}


PetscErrorCode IceHEINOModel::heinoCreateDat() {
  PetscErrorCode  ierr;
  char            tsname[2][4], tssname[3][5], tspname[3][4];
  char            filename[40],basefilename[30];

/* files are "PISMHEINO_##_??.dat" where ##=ismipRunName and ?? is one of the 
   following:
                  ts_iv,   ts_tba, 

                  tss_ait, tss_ahbt, tss_tba,

                  tsp1_it, tsp1_hbt, tsp1_bfh,
                  tsp2_it, tsp2_hbt, tsp2_bfh,
                  tsp3_it, tsp3_hbt, tsp3_bfh,
                  tsp4_it, tsp4_hbt, tsp4_bfh,
                  tsp5_it, tsp5_hbt, tsp5_bfh,
                  tsp6_it, tsp6_hbt, tsp6_bfh,
                  tsp7_it, tsp7_hbt, tsp7_bfh
*/

  strcpy(basefilename,heinodatprefix);
  strcat(basefilename,"_");
  strcat(basefilename,ismipRunName);
  ierr = PetscPrintf(grid.com,
            "creating ISMIP-HEINO deliverable files named %s_???.dat ...\n",
            basefilename); CHKERRQ(ierr);

  strcpy(tsname[0],"iv");
  strcpy(tsname[1],"tba");
  strcpy(tssname[0],"ait");
  strcpy(tssname[1],"ahbt");
  strcpy(tssname[2],"tba");
  strcpy(tspname[0],"it");
  strcpy(tspname[1],"hbt");
  strcpy(tspname[2],"bfh");
  
  // create global time series viewers (i.e. ASCII .dat files)
  for (PetscInt i=0; i<2; i++) {
    strcpy(filename,basefilename);
    strcat(filename,"_ts_");
    strcat(filename,tsname[i]);
    strcat(filename,".dat");
    ierr = vPetscPrintf(grid.com, "  opening %s\n", filename); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(grid.com, filename, &ts[i]); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(ts[i], PETSC_VIEWER_ASCII_DEFAULT); CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(ts[i],"SAMPLE JUNK in file %s\n", filename);  CHKERRQ(ierr);
  }

  // create sediment region time series viewers (i.e. ASCII .dat files)
  for (PetscInt i=0; i<3; i++) {
    strcpy(filename,basefilename);
    strcat(filename,"_tss_");
    strcat(filename,tssname[i]);
    strcat(filename,".dat");
    ierr = vPetscPrintf(grid.com, "  opening %s\n", filename); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(grid.com, filename, &tss[i]); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(tss[i], PETSC_VIEWER_ASCII_DEFAULT); CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(tss[i],"SAMPLE JUNK in file %s\n", filename);  CHKERRQ(ierr);
  }

  // create P1,...,P7 time series viewers (i.e. ASCII .dat files)
  for (PetscInt j=0; j<7; j++) {
    for (PetscInt i=0; i<3; i++) {
      strcpy(filename,basefilename);
      strcat(filename,"_tsp");
      const char nums[]="1234567";
      char       foo[]="A";
      foo[0] = nums[j];
      strcat(filename,foo);
      strcat(filename,"_");
      strcat(filename,tspname[i]);
      strcat(filename,".dat");
      ierr = vPetscPrintf(grid.com, "  opening %s\n", filename); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(grid.com, filename, &tsp[j][i]); CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(tsp[j][i], PETSC_VIEWER_ASCII_DEFAULT); CHKERRQ(ierr);
//      ierr = PetscViewerASCIIPrintf(tsp[j][i],"SAMPLE JUNK in file %s\n", filename);  CHKERRQ(ierr);
    }
  }

  return 0;
}


PetscErrorCode IceHEINOModel::heinoCloseDat() {
  PetscErrorCode  ierr;

  for (PetscInt i=0; i<2; i++) {
    ierr = PetscViewerDestroy(ts[i]); CHKERRQ(ierr);
  }
  for (PetscInt i=0; i<3; i++) {
    ierr = PetscViewerDestroy(tss[i]); CHKERRQ(ierr);
  }
  for (PetscInt j=0; j<7; j++) {
    for (PetscInt i=0; i<3; i++) {
      ierr = PetscViewerDestroy(tsp[j][i]); CHKERRQ(ierr);
    }
  }
  return 0;
}


bool IceHEINOModel::nearHeino(const PetscScalar x1, const PetscScalar y1,
                            const PetscScalar x2, const PetscScalar y2) {
  const PetscScalar epsSuff = 1.0; // within epsSuff meters
  return ( (fabs(x1 - x2) < epsSuff) && (fabs(y1 - y2) < epsSuff) );
}


PetscErrorCode IceHEINOModel::additionalAtStartTimestep() {
  // this is called at the beginning of time-stepping loop in IceModel::run()

  if (getExperName() == '0') {
    if ((ismipNoDeliver == PETSC_TRUE) && (ismipAllowAdapt == PETSC_TRUE)) {
      // do nothing with time step; allow long time step
    } else if ((ismipNoDeliver == PETSC_FALSE) && (ismipAllowAdapt == PETSC_TRUE)) {
      // go to next integer year
      maxdt_temporary = (1.0 - fmod(grid.p->year, 1.0)) * secpera;
    } else {
      dt_force = 0.25 * secpera;  // totally override adaptive time-stepping
                                  // and use .25 year steps as spec-ed by HEINO
    }
  }
  return 0;
}


PetscErrorCode IceHEINOModel::additionalAtEndTimestep() {
  PetscErrorCode  ierr;
  // this is called at the end of time-stepping loop in IceModel::run()
  // here we compute the ISMIP-HEINO "deliverables" at every time step

  // format specified by ISMIP-HEINO is following FORTRAN:
  //       write(file_number, 1000) time, quantity
  //       1000 format(2e14.6)
  // which I will take as
  //      ierr = PetscViewerASCIIPrintf(ts[i],"%14.6E%14.6E\n", time, quantity);
  //              CHKERRQ(ierr);
  // for now    

  if ((getExperName() == '0') && (ismipNoDeliver != PETSC_TRUE)) {

    // only deliver at integer year; must be on the year to 9 digits before we
    // will proceed to deliver
    if ( (fmod(grid.p->year,1.0) / PetscMax(1.0,grid.p->year)) > 5.0e-10)
      return 0;
    
    // get volume in 10^6 km^3 
    // and temperate (at pressure-melting = (T > DEFAULT_MIN_TEMP_FOR_SLIDING))
    // basal area in 10^6 km^2
    PetscScalar     **H, **ub, **vb, ***T;
    PetscScalar     volume=0.0, meltedarea=0.0,
                    ssavthick=0.0, ssavbasetemp=0.0, ssmeltedarea=0.0;
    PetscScalar     gvolume, gmeltedarea,
                    gssavthick, gssavbasetemp, gssmeltedarea;
    PetscScalar     pnthick[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                    pnbasetemp[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                    pnbaseheating[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                    gpnthick[7], gpnbasetemp[7], gpnbaseheating[7];
    // x-coords (in ISMIP-HEINO scheme) for delivered grid points; note Py=2000e3 for all
    const PetscScalar  Px[7] = {3900e3, 3800e3, 3700e3, 3500e3, 3200e3, 2900e3, 2600e3};
    const PetscScalar  a = grid.p->dx * grid.p->dy * 1e-3 * 1e-3; // area unit (km^2)
  
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar dv = a * H[i][j] * 1e-3; // column vol (km^3)
        volume += dv;
        const PetscScalar homT0 = T[i][j][0] + ice.beta_CC_grad * H[i][j];
        const bool ismelted = ((H[i][j] > 0) && (homT0 > DEFAULT_MIN_TEMP_FOR_SLIDING));
        if (ismelted)
          meltedarea += a;
        const PetscScalar xIH = grid.p->dx * i, 
                          yIH = grid.p->dy * j;
        if (inSoftSediment(xIH,yIH)) {
          ssavthick += H[i][j] * a;  // see averaging below to remove area factor
          ssavbasetemp += homT0 * a;
          if (ismelted)
            ssmeltedarea += a;
          for (PetscInt k=0; k<7; k++) {
            if (nearHeino(xIH,yIH,Px[k],2000e3)) {
              // ierr = PetscPrintf(grid.com,
              //          "  [at ISMIP-HEINO point P%d]\n",k+1); CHKERRQ(ierr);
              pnthick[k] = H[i][j];
              pnbasetemp[k] = homT0;
              // note equation (10) in ISMIP-HEINO documentation [date: 19 July 2006]
              // is wrong; should read
              //     " \tau_b = - \rho g H \grad_H h " 
              const PetscScalar 
                  taubx = -ice.rho * ice.grav * H[i][j]
                            * (H[i+1][j] - H[i-1][j]) / (2 * grid.p->dx),
                  tauby = -ice.rho * ice.grav * H[i][j] 
                            * (H[i][j+1] - H[i][j-1]) / (2 * grid.p->dy);
              pnbaseheating[k] = taubx * ub[i][j] + tauby * vb[i][j];
//              if (pnbaseheating[k] < 0.0) {
//                SETERRQ(1, 
//                  "basal heating negative in IceHEINOModel::additionalStuffAtTimestep()\n");
//              }
            }
          }
        }
      }
    }  
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

    ierr = PetscGlobalSum(&volume, &gvolume, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&meltedarea, &gmeltedarea, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&ssavthick, &gssavthick, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&ssavbasetemp, &gssavbasetemp, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&ssmeltedarea, &gssmeltedarea, grid.com); CHKERRQ(ierr);
    // actually average
    PetscScalar ssTotalArea = 1000e3 + 140e3; // 1,140,000 (km^2)
    const PetscScalar dxkm = grid.p->dx / 1000.0, dykm = grid.p->dy / 1000.0;
    ssTotalArea += 2 * 1700.0 * (dykm/2) + 2 * 1000.0 * (dxkm/2)
                   + 4 * (dxkm/2) * (dykm/2); // correct for edge effect;
    gssavthick = gssavthick / ssTotalArea;
    gssavbasetemp = gssavbasetemp / ssTotalArea;

    ierr = PetscViewerASCIIPrintf(ts[0],"%14.6E%14.6E\n", 
               grid.p->year, gvolume * 1e-6); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(ts[1],"%14.6E%14.6E\n", 
               grid.p->year, gmeltedarea * 1e-6); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(tss[0],"%14.6E%14.6E\n", 
               grid.p->year, gssavthick * 1e-3); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(tss[1],"%14.6E%14.6E\n", 
               grid.p->year, gssavbasetemp); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(tss[2],"%14.6E%14.6E\n", 
               grid.p->year, gssmeltedarea * 1e-6); CHKERRQ(ierr);

    for (PetscInt k=0; k<7; k++) {
      ierr = PetscGlobalMax(&pnthick[k], &gpnthick[k], grid.com); CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(tsp[k][0],"%14.6E%14.6E\n", 
               grid.p->year, gpnthick[k] * 1e-3); CHKERRQ(ierr);
      ierr = PetscGlobalMax(&pnbasetemp[k], &gpnbasetemp[k], grid.com); CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(tsp[k][1],"%14.6E%14.6E\n", 
               grid.p->year, gpnbasetemp[k]); CHKERRQ(ierr);
      ierr = PetscGlobalMax(&pnbaseheating[k], &gpnbaseheating[k], grid.com); CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(tsp[k][2],"%14.6E%14.6E\n", 
               grid.p->year, gpnbaseheating[k]); CHKERRQ(ierr);
    }

  }
  return 0;
}


PetscErrorCode IceHEINOModel::simpFinalize() {
  PetscErrorCode  ierr;
  if (ismipNoDeliver != PETSC_TRUE) {
    ierr = heinoCloseDat(); CHKERRQ(ierr);
  }
  return 0;
}
