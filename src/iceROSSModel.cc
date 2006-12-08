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
  ierr = VecDuplicate(vh, &ubarBC); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbarBC); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::destroyROSSVecs() {
  PetscErrorCode ierr;

  ierr = VecDestroy(obsAzimuth); CHKERRQ(ierr);
  ierr = VecDestroy(obsMagnitude); CHKERRQ(ierr);
  ierr = VecDestroy(obsAccurate); CHKERRQ(ierr);
  ierr = VecDestroy(ubarBC); CHKERRQ(ierr);
  ierr = VecDestroy(vbarBC); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::initFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      ROSSchosen, inFileSet, prefixSet;
  char            inFile[PETSC_MAX_PATH_LEN];

//   This option chooses EISMINT ROSS, i.e. from the paper
//     MacAyeal and five others (1996). "An ice-shelf model test based on the Ross ice shelf,"
//     Ann. Glaciol. 23, 46-51
  ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &ROSSchosen); CHKERRQ(ierr);
  if (ROSSchosen == PETSC_FALSE) {
    SETERRQ(1,"why did I get here? (IceROSSModel::initFromOptions)");
  }
  
  ierr = verbPrintf(2,grid.com, 
            "initializing EISMINT ROSS validation test ... \n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    SETERRQ(2,"PISM input file not allowed for initialization of EISMINT ROSS");
  }
  
  ierr = PetscOptionsGetString(PETSC_NULL, "-prefix", prefixROSS,
                               PETSC_MAX_PATH_LEN, &prefixSet); CHKERRQ(ierr);
  if (prefixSet == PETSC_FALSE) {
    strcpy(prefixROSS,"eisROSS/");
  }
  ierr = verbPrintf((prefixSet == PETSC_TRUE) ? 2 : 4, grid.com,
            "prefix for ROSS data files =%s\n",prefixROSS); CHKERRQ(ierr);
  
  ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);

  // note region is *roughly* 1000km by 750km, but MacAyeal et al 1996 specifies these:
  PetscScalar dxROSS = 6822;
  // expand grid in y direction, but also swap meaning of x and y
  xsROSS = 36; xmROSS = 111;
  PetscInt    MxROSS = xsROSS + xmROSS, MyROSS = 147;  // both are 147; square grid
  grid.p->Mx=MxROSS;
  grid.p->My=MyROSS;
  // Mz is 31 by default but can be set by user
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = grid.rescale(0.5 * (MxROSS - 1) * dxROSS, 0.5 * (MyROSS - 1) * dxROSS, 1000); CHKERRQ(ierr);

  ierr = createROSSVecs(); CHKERRQ(ierr);

  ierr = readROSSfiles(); CHKERRQ(ierr); // reads from 111by147Grid.dat and kbc.dat
  
  ierr = fillinTemps();  CHKERRQ(ierr);
  
  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running EISMINT ROSS ...\n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceROSSModel::fillinTemps() {
  PetscErrorCode      ierr;
  PetscScalar         **Ts, ***T, ***Tb, ***tau;

  // fill in all temps with Ts
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<grid.p->Mz; k++) {
        T[i][j][k] = Ts[i][j];
        tau[i][j][k] = 0.0;
      }
      for (PetscInt k=0; k<grid.p->Mbz; k++)
        Tb[i][j][k] = Ts[i][j];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::readROSSfiles() {
  // reads from 
  //               111by147Grid.dat
  //               kbc.dat
  //               inlets.dat
  PetscErrorCode  ierr;
  char            datfilename[PETSC_MAX_PATH_LEN], kbcfilename[PETSC_MAX_PATH_LEN],
                  inletsfilename[PETSC_MAX_PATH_LEN];
  char            ignor[400];
  PetscScalar     **mask, **azimuth, **magnitude, **H, 
                  **accurate, **bed, **accum, **Ts;
  
  // some grid info
  ierr = verbPrintf(5,grid.com, 
            "  readROSSFiles(): grid.xs = %d, grid.xm = %d, grid.ys = %d, grid.ym = %d\n",
            grid.xs, grid.xm, grid.ys, grid.ym); CHKERRQ(ierr);
  ierr = verbPrintf(5,grid.com, 
            "  readROSSFiles(): xsROSS = %d, xmROSS = %d\n",xsROSS,xmROSS); CHKERRQ(ierr);

  strcpy(datfilename,prefixROSS);
  strcat(datfilename,"111by147Grid.dat");
  strcpy(kbcfilename,prefixROSS);
  strcat(kbcfilename,"kbc.dat");
  strcpy(inletsfilename,prefixROSS);
  strcat(inletsfilename,"inlets.dat");
  
  ierr = verbPrintf(2,grid.com,"reading from data file %s ...\n",datfilename); CHKERRQ(ierr);
  FILE* datfile = fopen(datfilename,"r");
  if (datfile == NULL) {
    SETERRQ1(1,"error opening file %s for reading in IceROSSModel",datfilename);
  }
  
  // ignor first 4 lines, which contain comments
  for (PetscInt k = 0; k < 4; k++) {
    fgets(ignor, 400, datfile);
    if (ferror(datfile)) {
      SETERRQ1(1,"error reading from file %s in IceROSSModel",datfilename);
    }
  }

  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines give latitude (WHAT SYSTEM?) of rows
  for (PetscInt i=0; i < xmROSS; i++) {
    float lat;
    fscanf(datfile, "%f", &lat);
    gridLat[i] = (PetscScalar) lat;  // latitude of row H[i+xsROSS][:], for ex.
  }

  // ignor another couple of lines
  fgets(ignor, 400, datfile);  // actually just gets end of line ...
  fgets(ignor, 400, datfile);  //HUH?  WHY NECESSARY?
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 147 lines give longitude (WHAT SYSTEM?) of rows
  for (PetscInt j=0; j < 147; j++) {
    float lon;
    fscanf(datfile, "%f", &lon);
    gridLon[j] = (PetscScalar) lon;  // longitude of column H[:][j], for ex.
  }

  // ignor another couple of lines
  fgets(ignor, 400, datfile);  // actually just gets end of line ...
  fgets(ignor, 400, datfile);  //HUH?  WHY NECESSARY?
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives mask (initially; region B not right)
  ierr = VecSet(vMask, MASK_FLOATING); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      int mm;
      fscanf(datfile, "%d", &mm);
      mask[i][j] = (mm == 1) ? MASK_FLOATING : MASK_SHEET;
    }
  }
  // write in two rows of MASK_SHEET to "block" interaction of top and bottom
  for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
    mask[0][j] = MASK_SHEET;
    mask[1][j] = MASK_SHEET;
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  // ignor another couple of lines
  fgets(ignor, 400, datfile);  // actually just gets end of line ...
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);
 
  // next 111 lines gives azimuth (of velocity)
  ierr = VecSet(obsAzimuth, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azimuth); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      float azi;
      fscanf(datfile, "%f", &azi);
      azimuth[i][j] = (PetscScalar) azi;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, obsAzimuth, &azimuth); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives magnitude (of velocity); also fill in velocity
  ierr = VecSet(obsMagnitude, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &magnitude); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      float mag;
      fscanf(datfile, "%f", &mag);
      magnitude[i][j] = ((PetscScalar) mag) / secpera;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, obsMagnitude, &magnitude); CHKERRQ(ierr);    

  // ignor another couple of lines
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  fgets(ignor, 400, datfile);
  ierr = verbPrintf(4,grid.com,"  LAST IGNORED LINE READS:%s", ignor); CHKERRQ(ierr);

  // next 111 lines gives thickness
  ierr = VecSet(vH, 1.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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
  ierr = VecSet(obsAccurate, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, obsAccurate, &accurate); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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
  ierr = VecSet(vbed, -600.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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
  ierr = VecSet(vAccum, 0.2 / secpera); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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
  ierr = VecSet(vTs, 248.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
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

  /*****************************************************************/
  /*       NOW READ   kbc.dat   and use to set b.c.s for velocity  */
  /*****************************************************************/

  PetscScalar   **ubar, **vbar;  // meaning switched!
  ierr = VecSet(ubarBC,0.0); CHKERRQ(ierr);
  ierr = VecSet(vbarBC,0.0); CHKERRQ(ierr);  
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, ubarBC, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vbarBC, &vbar); CHKERRQ(ierr);    

  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azimuth); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &magnitude); CHKERRQ(ierr);    
  ierr = verbPrintf(2,grid.com,"reading from data file %s ...\n",kbcfilename); CHKERRQ(ierr);
  FILE* kbcfile = fopen(kbcfilename,"r");
  if (kbcfile == NULL) {
    SETERRQ1(1,"error opening file %s for reading in IceROSSModel",kbcfilename);
  }
  const PetscScalar pi = 3.14159265358979;
  for (PetscInt k=0; k < 77; k++) {
    int temp_i, temp_j;
    fscanf(kbcfile, "%d %d", &temp_i, &temp_j);
    const PetscInt i = ((PetscInt) temp_i)+xsROSS, 
                   j = (PetscInt) temp_j;
    kbcGridLoc[0][k] = i;
    kbcGridLoc[1][k] = j;
    vbar[i][j] = magnitude[i][j] * sin((pi/180.0) * azimuth[i][j]);
    ubar[i][j] = magnitude[i][j] * cos((pi/180.0) * azimuth[i][j]);
    mask[i][j] = MASK_SHEET;  // so that they are fixed
  }
  if (fclose(kbcfile) != 0) {
    SETERRQ1(1,"error closing file %s in IceROSSModel",kbcfilename);
  }
  ierr = DAVecRestoreArray(grid.da2, obsAzimuth, &azimuth); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, obsMagnitude, &magnitude); CHKERRQ(ierr);    

  ierr = verbPrintf(2,grid.com,"reading from data file %s ...\n",inletsfilename); CHKERRQ(ierr);
  FILE* inlfile = fopen(inletsfilename,"r");
  if (inlfile == NULL) {
    SETERRQ1(1,"error opening file %s for reading in IceROSSModel",inletsfilename);
  }
  for (PetscInt k=0; k < 22; k++) {
    int            temp_i, temp_j;
    float          temp_azi, temp_mag;
    fscanf(inlfile, "%d %d %f %f", &temp_i, &temp_j, &temp_azi, &temp_mag);
    const PetscInt i = ((PetscInt) temp_i)+xsROSS, 
                   j = (PetscInt) temp_j;
    const PetscScalar mag = ((PetscScalar) temp_mag) / secpera;
    inletGridLoc[0][k] = i;
    inletGridLoc[1][k] = j;
    vbar[i][j] = mag * sin((pi/180.0) * ((PetscScalar) temp_azi));
    ubar[i][j] = mag * cos((pi/180.0) * ((PetscScalar) temp_azi));
    mask[i][j] = MASK_SHEET;  // so that they are fixed
  }
  if (fclose(inlfile) != 0) {
    SETERRQ1(1,"error closing file %s in IceROSSModel",inletsfilename);
  }

  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, ubarBC, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vbarBC, &vbar); CHKERRQ(ierr);    
  return 0;
}


PetscErrorCode IceROSSModel::setBoundaryVels() {
  PetscErrorCode ierr;
  PetscScalar   **ubar, **vbar, **uvbar[2];  // meaning switched!
  
  // set initial velocities in shelf for iteration
  ierr = VecSet(vubar,-200.0 / secpera); CHKERRQ(ierr);
  ierr = VecSet(vvbar,0.0); CHKERRQ(ierr);

  // set boundary condition which will apply to finite difference system:
  //    staggered grid velocities at MASK_SHEET points which neighbor MASK_FLOATING points
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, ubarBC, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vbarBC, &vbar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  for (PetscInt k=0; k < 77; k++) {
    const PetscInt i = kbcGridLoc[0][k];
    const PetscInt j = kbcGridLoc[1][k];
    uvbar[1][i][j-1] = vbar[i][j];
    uvbar[1][i][j] = vbar[i][j];
    uvbar[0][i-1][j] = ubar[i][j];
    uvbar[0][i][j] = ubar[i][j];
  }
  for (PetscInt k=0; k < 22; k++) {
    const PetscInt i = inletGridLoc[0][k];
    const PetscInt j = inletGridLoc[1][k];
    uvbar[1][i][j-1] = vbar[i][j];
    uvbar[1][i][j] = vbar[i][j];
    uvbar[0][i-1][j] = ubar[i][j];
    uvbar[0][i][j] = ubar[i][j];
  }
  ierr = DAVecRestoreArray(grid.da2, ubarBC, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vbarBC, &vbar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceROSSModel::makeSurfaceFloating() {
  PetscErrorCode  ierr;
  PetscScalar     **h, **H;
  // update surface elev: APPLY FLOATING CRITERION EVERYWHERE to get smooth surface
  ierr = VecSet(vh, (1-ice.rho/ocean.rho) * 1.0); CHKERRQ(ierr);  // start with all 1m case
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      h[i][j] = (1-ice.rho/ocean.rho) * H[i][j];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);    

  return 0;
}


PetscErrorCode IceROSSModel::showObservedVels() {
  PetscErrorCode  ierr;
  PetscScalar **azi, **mag, **acc, **ubar, **vbar;   
  const PetscScalar pi = 3.14159265358979;
  ierr = verbPrintf(2,grid.com,
              "updating viewers to show observed velocities ...\n"); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azi); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &mag); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsAccurate, &acc); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      vbar[i][j] = acc[i][j] * mag[i][j] * sin((pi/180.0) * azi[i][j]);
      ubar[i][j] = acc[i][j] * mag[i][j] * cos((pi/180.0) * azi[i][j]);
      // uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
      // uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))
    }
  }
  ierr = DAVecRestoreArray(grid.da2, obsMagnitude, &mag); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, obsAzimuth, &azi); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, obsAccurate, &acc); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);

  ierr = updateViewers(); CHKERRQ(ierr);

  Vec vNu[2] = {vWork2d[0], vWork2d[1]};
  ierr = computeEffectiveViscosity(vNu, macayealEpsilon); CHKERRQ(ierr);
  ierr = updateNuViewers(vNu,vNu,false); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::writeROSSfiles() {
  PetscErrorCode  ierr;
  PetscTruth      outputSet, formatSet;
  char            basename[PETSC_MAX_PATH_LEN] = "pismROSS", format[10];
  // writes h, H, bed, mask, and c=sqrt(ubar^2 + vbar^2) into pismROSS.m
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", basename,
                               PETSC_MAX_PATH_LEN, &outputSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-of", format,
                               10, &formatSet); CHKERRQ(ierr);
  if ((outputSet == PETSC_TRUE) && (formatSet == PETSC_TRUE)
      && (strchr(format,'m') != NULL)) {
    ierr = writeFiles(basename,"m"); CHKERRQ(ierr);  
    ierr = verbPrintf(2,grid.com, "\n"); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceROSSModel::readRIGGSandCompare() {
  PetscErrorCode  ierr;
  char            riggsfilename[PETSC_MAX_PATH_LEN];
  char            ignor[400];

  strcpy(riggsfilename,prefixROSS);
  strcat(riggsfilename,"RIGGS.dat");
  ierr = verbPrintf(2,grid.com,"reading from data file %s ...\n",riggsfilename); CHKERRQ(ierr);
  FILE* riggsfile = fopen(riggsfilename,"r");
  if (riggsfile == NULL) {
    SETERRQ1(1,"error opening file %s for reading in IceROSSModel",riggsfilename);
  }
  // ignor first 4 lines;  JUST FOR TEST OF READ FOR NOW
  for (PetscInt j = 1; j <= 4; j++) {
    fgets(ignor, 400, riggsfile);
    if (ferror(riggsfile)) {
      SETERRQ1(1,"error reading from file %s in IceROSSModel",riggsfilename);
    }
  }
  if (fclose(riggsfile) != 0) {
    SETERRQ1(1,"error closing file %s in IceROSSModel",riggsfilename);
  }
  return 0;
}


PetscErrorCode IceROSSModel::computeErrorsInAccurate() {
  // average over grid, where observed velocity is accurate according to 111by147grid.dat,
  // the difference between computed and observed u,v
  PetscErrorCode  ierr;
  PetscScalar  uerr=0.0, verr=0.0,
               accN=0.0, accArea=0.0, maxcComputed=0.0,
               **azi, **mag, **acc, **ubar, **vbar, **H, **mask;
  vecErrAcc = 0.0;
  const PetscScalar pi = 3.14159265358979, area = grid.p->dx * grid.p->dy;
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azi); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &mag); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsAccurate, &acc); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  for (PetscInt i=xsROSS; i < xsROSS+xmROSS; i++) {
    for (PetscInt j=grid.ys; j < grid.ys+grid.ym; j++) {
      if ((modMask(mask[i][j]) == MASK_FLOATING) && (H[i][j] > 1.0)) {
        const PetscScalar ccomputed = sqrt(PetscSqr(vbar[i][j]) + PetscSqr(ubar[i][j]));
        maxcComputed = PetscMax(maxcComputed,ccomputed);
        if (PetscAbs(acc[i][j] - 1.0) < 0.1) {
          accN += 1.0;
          accArea += area;
          const PetscScalar vobs = mag[i][j] * sin((pi/180.0) * azi[i][j]);
          const PetscScalar uobs = mag[i][j] * cos((pi/180.0) * azi[i][j]);
          // compare from readme.txt:
          // uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
          // uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))
          const PetscScalar Dv = PetscAbs(vobs - vbar[i][j]);
          const PetscScalar Du = PetscAbs(uobs - ubar[i][j]);
          verr += Dv;
          uerr += Du;
          vecErrAcc += (PetscSqr(Dv) + PetscSqr(Du)) * area;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, obsMagnitude, &mag); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, obsAzimuth, &azi); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, obsAccurate, &acc); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  // are sums of m/s; report in m/a and as averages
  ierr = verbPrintf(2,grid.com,"maximum computed speed in ice shelf is %10.3f (m/a)\n",
             maxcComputed * secpera); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"ERRORS relative to observations of Ross Ice Shelf:\n"); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  [number of grid points in, and area of, 'accurate observed area' are %d, %9.4f (km^2)]\n",
             (int) accN, accArea / 1e6); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average (over accurate area) error in x-comp of vel   = %9.3f (m/a)\n",
             (verr * secpera) / accN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average (over accurate area) error in y-comp of vel   = %9.3f (m/a)\n",
             (uerr * secpera) / accN); CHKERRQ(ierr);
  vecErrAcc = secpera * sqrt(vecErrAcc) / sqrt(accArea);
  ierr = verbPrintf(2,grid.com,
             "  rms average (over accurate area) error in vector vel  = %9.3f (m/a)\n",
             vecErrAcc); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::runTune() {
  PetscErrorCode  ierr;
  PetscReal       tuneparam[20];
  PetscTruth      tuneSet;
  PetscInt        numtuneparam=3;
  
  // user might enter
  //     obj/pisms -ross -tune 1.7e8,0.1e8,2.5e8
  // [note no spaces between parameters!]
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-tune", tuneparam, &numtuneparam, &tuneSet); CHKERRQ(ierr);
  if (tuneSet == PETSC_FALSE || numtuneparam != 3) {
    SETERRQ(1,"-tune not accompanied by three parameters");
  }
  
  ierr = makeSurfaceFloating(); CHKERRQ(ierr);

  useMacayealVelocity = PETSC_TRUE;
  useConstantNuForMacAyeal = PETSC_FALSE;
  useConstantHardnessForMacAyeal = PETSC_TRUE;
  setMacayealEpsilon(0.0);  // don't use this lower bound

  regularizingVelocitySchoof = 1.0 / secpera;  // (VELOCITY/LENGTH)^2  is very close to 10^-27
  regularizingLengthSchoof = 1000.0e3;
  ierr = verbPrintf(2,grid.com,"[using Schoof regularization constant = %10.5e for all runs]\n",
              PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof)); CHKERRQ(ierr);

  // in Matlab notation:  hard = tuneparam[0]:tuneparam[1]:tuneparam[2]
  for (PetscScalar hard = tuneparam[0]; hard < tuneparam[2]+0.1*tuneparam[1]; hard += tuneparam[1]) { 
    ierr = verbPrintf(2,grid.com,"computing velocities with hardness (=bar B) = %10.5e ...\n",hard); CHKERRQ(ierr);
    ierr = setBoundaryVels(); CHKERRQ(ierr);
    constantHardnessForMacAyeal = hard;
    // solve ice shelf equations:
    ierr = velocityMacayeal(); CHKERRQ(ierr);
    ierr = updateViewers(); CHKERRQ(ierr);
    // compare to observed
    ierr = computeErrorsInAccurate(); CHKERRQ(ierr);
  }

  ierr = writeROSSfiles(); CHKERRQ(ierr);
  ierr = destroyROSSVecs(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::run() {
  PetscErrorCode  ierr;
  PetscInt        pause_time;
  PetscTruth      pause_p, showobsvel, tune;

  ierr = PetscOptionsHasName(PETSC_NULL, "-tune", &tune); CHKERRQ(ierr);
  if (tune == PETSC_TRUE) {
    ierr = runTune(); CHKERRQ(ierr);
    return 0;
  }

  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, &pause_p); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-showobsvel", &showobsvel); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
  "$$$$$      YEAR (+    STEP[$]):     VOL    AREA    MELTF     THICK0     TEMP0\n");
  CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep!

  ierr = makeSurfaceFloating(); CHKERRQ(ierr);

  ierr = setBoundaryVels(); CHKERRQ(ierr);

  // solve model equations 
  useMacayealVelocity = PETSC_TRUE;
  useConstantNuForMacAyeal = PETSC_FALSE;
  useConstantHardnessForMacAyeal = PETSC_TRUE;
  setMacayealEpsilon(0.0);  // don't use this lower bound
//  constantHardnessForMacAyeal = 1.9e8;  // Pa s^{1/3}
  constantHardnessForMacAyeal = 2.22e8;  // Pa s^{1/3}

  regularizingVelocitySchoof = 1.0 / secpera;  // (VELOCITY/LENGTH)^2  is very close to 10^-27
  regularizingLengthSchoof = 1000.0e3;
  ierr = verbPrintf(5,grid.com,"[using Schoof regularization constant = %10.5e]\n",
              PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof)); CHKERRQ(ierr);

  ierr = velocityMacayeal(); CHKERRQ(ierr);

  // report on result of computation (i.e. to standard out, to viewers, to Matlab file)
  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  ierr = summary(true,true); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);
  // compare to observed
  ierr = computeErrorsInAccurate(); CHKERRQ(ierr);
  // ierr = readRIGGSandCompare(); CHKERRQ(ierr);
  ierr = writeROSSfiles(); CHKERRQ(ierr);
  if (pause_p == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }

  if (showobsvel == PETSC_TRUE) { 
    // update viewers to show observed velocities
    ierr = showObservedVels(); CHKERRQ(ierr);
    if (pause_p == PETSC_TRUE) {
      ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
      ierr = PetscSleep(pause_time); CHKERRQ(ierr);
    }
  }

  ierr = destroyROSSVecs(); CHKERRQ(ierr);
  return 0;
}
