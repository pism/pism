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
#include <cstdio>
#include <cmath>
#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

#include "iceROSSModel.hh"


IceROSSModel::IceROSSModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {  // do nothing; note derived classes must have constructors
}


void IceROSSModel::setflowlawNumber(PetscInt law) {
  flowlawNumber = law;
}


PetscInt IceROSSModel::getflowlawNumber() {
  return flowlawNumber;
}


PetscErrorCode IceROSSModel::createROSSVecs() {
  PetscErrorCode ierr;

  ierr = VecDuplicate(vh, &obsAzimuth); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &obsMagnitude); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &obsAccurate); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::destroyROSSVecs() {
  PetscErrorCode ierr;

  ierr = VecDestroy(obsAzimuth); CHKERRQ(ierr);
  ierr = VecDestroy(obsMagnitude); CHKERRQ(ierr);
  ierr = VecDestroy(obsAccurate); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::initFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      ROSSchosen, inFileSet;
  char            inFile[PETSC_MAX_PATH_LEN];

//   This option chooses EISMINT ROSS, i.e. from the paper
//     MacAyeal and five others (1996). "An ice-shelf model test based on the Ross ice shelf,"
//     Ann. Glaciol. 23, 46-51
  ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &ROSSchosen); CHKERRQ(ierr);
  if (ROSSchosen == PETSC_FALSE) {
    SETERRQ(1,"why did I get here? (IceROSSModel::initFromOptions)");
  }
  
  ierr = PetscPrintf(grid.com, 
            "initializing EISMINT ROSS validation test ... \n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    SETERRQ(2,"PISM input file not allowed for initialization of EISMINT ROSS");
  }
  
  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);

  // note region is *roughly* 1000km by 750km, but MacAyeal et al 1996 specifies these:
  PetscScalar dxROSS = 6822;
  PetscInt    MxROSS = 111, MyROSS = 147;  // swap from usual
  grid.p->Mx=MxROSS;
  grid.p->My=MyROSS;
  // Mz is 31 by default but can be set by user
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = grid.rescale(0.5 * (MxROSS - 1) * dxROSS, 0.5 * (MyROSS - 1) * dxROSS, 1000); CHKERRQ(ierr);

  ierr = createROSSVecs(); CHKERRQ(ierr);

  ierr = readROSSfile(); CHKERRQ(ierr); // reads from 111by147Grid.dat  !!
  
  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = PetscPrintf(grid.com, "running EISMINT ROSS ...\n"); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode IceROSSModel::readROSSfile() {
  // reads from eisROSS/111by147Grid.dat  !!
  PetscErrorCode  ierr;
  char            datfilename[PETSC_MAX_PATH_LEN] = "eisROSS/111by147Grid.dat";
  char            ignor[400];
  const PetscScalar pi = 3.14159265358979;
  PetscScalar     **mask, **azimuth, **magnitude, **ubar, **vbar, **H, 
                  **accurate, **bed, **accum, **Ts;

  FILE* datfile = fopen(datfilename,"r");
  
  if (datfile == NULL) {
    SETERRQ1(1,"error opening file %s for reading in IceROSSModel",datfilename);
  }
  
  // ignor first 268 lines, which contain comments and lat/long of locations
  for (PetscInt j = 1; j <= 268; j++) {
    fgets(ignor, 400, datfile);
    if (ferror(datfile)) {
      SETERRQ1(1,"error reading from file %s in IceROSSModel",datfilename);
    }
  }

//  ierr = verbPrintf(1,grid.com, 
//            "  readROSSFiles(): grid.xs = %d, grid.xm = %d, grid.ys = %d, grid.ym = %d\n",
//            grid.xs, grid.xm, grid.ys, grid.ym); CHKERRQ(ierr);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives mask (initially; region B not right)
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      int mm;
      fscanf(datfile, "%d", &mm);
      mask[i][j] = (mm == 1) ? MASK_FLOATING : MASK_SHEET;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  // ignor another couple of lines
  fgets(ignor, 400, datfile);  // actually just gets end of line ...
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);
 
  // next 111 lines gives azimuth (of velocity)
  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azimuth); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      float azi;
      fscanf(datfile, "%f", &azi);
      azimuth[i][j] = (PetscScalar) azi;
    }
  }

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives magnitude (of velocity); also fill in velocity
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &magnitude); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      float mag;
      fscanf(datfile, "%f", &mag);
      magnitude[i][j] = ((PetscScalar) mag) / secpera;
      ubar[i][j] = magnitude[i][j] * sin((pi/180.0) * azimuth[i][j]);
      vbar[i][j] = magnitude[i][j] * cos((pi/180.0) * azimuth[i][j]);
      // uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
      // uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))
    }
  }
  ierr = DAVecRestoreArray(grid.da2, obsMagnitude, &magnitude); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, obsAzimuth, &azimuth); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives thickness
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float thick;
       fscanf(datfile, "%f", &thick);
       H[i][j] = (PetscScalar) thick;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives whether interpolation is accurate
  ierr = DAVecGetArray(grid.da2, obsAccurate, &accurate); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float acc;
       fscanf(datfile, "%f", &acc);
       accurate[i][j] = (PetscScalar) acc;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, obsAccurate, &accurate); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives sea bed depth; read it though not planned for use
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float depth;
       fscanf(datfile, "%f", &depth);
       bed[i][j] = - ((PetscScalar) depth);  // elevation is negative
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines determines whether to set thickness to 1 m
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float regionB;
       fscanf(datfile, "%f", &regionB);
       if (fabs(regionB - 1.0) < 0.1) {
         H[i][j] = 1.0;
       }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives accumulation in mm/yr
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float accum_mm;
       fscanf(datfile, "%f", &accum_mm);
       accum[i][j] = (((PetscScalar) accum_mm) / 1000.0) / secpera;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives \bar B; see comments in MacAyeal et al 1996 and in readme.txt
  // ignored for now
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float barB;
       fscanf(datfile, "%f", &barB);
    }
  }

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives surface temp in deg C
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i < grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
       float temp_C;
       fscanf(datfile, "%f", &temp_C);
       Ts[i][j] = ((PetscScalar) temp_C) + 273.15;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);    

  if (fclose(datfile) != 0) {
    SETERRQ1(1,"error closing file %s in IceROSSModel",datfilename);
  }

  return 0;
}



PetscErrorCode IceROSSModel::run() {
  PetscErrorCode  ierr;

  ierr = initSounding(); CHKERRQ(ierr);

  ierr = PetscPrintf(grid.com,
  "$$$$$      YEAR (+    STEP[$]):     VOL    AREA    MELTF     THICK0     TEMP0\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com, "$$$$$"); CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  ierr = summary(true,true); CHKERRQ(ierr);  // report starting state

  useMacayealVelocity = PETSC_TRUE;
  // ierr = velocity(true); CHKERRQ(ierr);
   
  ierr = summary(true,true); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);

  PetscInt pause_time;
  PetscTruth pause_p;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, &pause_p); CHKERRQ(ierr);
  if (pause_p == PETSC_TRUE)
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);

  ierr = destroyROSSVecs(); CHKERRQ(ierr);

  return 0;
}
