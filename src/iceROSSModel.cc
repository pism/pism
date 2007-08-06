// Copyright (C) 2006-2007 Ed Bueler
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
#include <cstdio>
#include <cmath>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

#include "iceROSSModel.hh"

/*  THIS COMMENT IS NO LONGER ACCURATE but it should be kept until new setup
is fully documented in the manual:
Run this derived class with "obj/pisms -ross".

Allows the following options:  
  -d cnmu       most useful way to see what is going on
  -if foo       NOT allowed!
  -o foo -of m  writes data to foo.m; note .pb is probably not useful
  -pause N      pause for N seconds when refreshing viewers
  -prefix foo   looks for files "111by147Grid.dat" and "kbc.dat" and 
                "inlets.dat" in foo/; files at
                http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html
  -ross         to start this
  -showobsvel   shows observed (but interpolated) speeds from "111by147Grid.dat"
  -tune x,y,z   run through \bar B=x:y:z (Matlab) or 
                \bar B = x, x+y, x+2y, ..., x+Ny=z as
                hardness parameters; note no spaces in "x,y,z"!
  -verbose      shows, in particular, omitted lines in "???.dat" reads

Example 1.  Basic run with all info displayed.  Note \bar B = 2.22e8 versus 
/bar B = 1.9e8 as in MacAyeal et al 1996.  Note files eisROSS/111by147Grid.dat,
eisROSS/kbc.dat, and eisROSS/inlets.dat must be present: 
  user@home:~/pism$ obj/pisms -ross -d cnmu -pause 10 -showobsvel -verbose
  
Example 2.  Same as above but saving Matlab results; no display.  The resulting
file can be viewed in Matlab.  See also test/ross/README.rossPISM, 
test/ross/riggs_ELBclean.dat, and test/ross/rossspeedplot.m, which can be used
to produce \Chi^2 statistic relative to RIGGS data, and nice picture: 
  user@home:~/pism$ obj/pisms -ross -verbose -o PISM_ross_2p22e8 -of m
  
Example 3.  Same, but asking for a lot more accuracy.  Nearly the same result, 
which suggests defaults (-mv_rtol 1e-4 -ksp_rtol 1e-6) suffice:
  user@home:~/pism$ obj/pisms -ross -d cnmu -pause 10 -showobsvel -verbose\
                    -mv_rtol 1e-7 -ksp_rtol 1e-10  

Example 4:  Tune across range of values of \bar B, including MacAyeal et al 1996
value and (close to) optimal value:
  user@home:~/pism$ obj/pisms -ross -verbose -tune 1.7e8,1e7,2.4e8
  
*/

IceROSSModel::IceROSSModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {  // do nothing; note derived classes must have constructors

  // use defined velocity boundary condition for shelf, instead of computing 
  // SIA values at those locations
  computeSIAVelocities = PETSC_FALSE;

  // further settings for velocity computation 
  useConstantNuForMacAyeal = PETSC_FALSE;
  useConstantHardnessForMacAyeal = PETSC_TRUE;
  setMacayealEpsilon(0.0);  // don't use this lower bound
  constantHardnessForMacAyeal = 1.9e8;  // Pa s^{1/3}; (MacAyeal et al 1996) value
  regularizingVelocitySchoof = 1.0 / secpera;  // 1 m/a is small velocity for shelf!
  regularizingLengthSchoof = 1000.0e3;         // (VELOCITY/LENGTH)^2  is very close to 10^-27
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

  ierr = IceModel::initFromOptions(PETSC_FALSE); CHKERRQ(ierr);

  // set Lz to 1000.0 m by default; note max thickness in data is 800.0 m
  PetscTruth LzSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-Lz", &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_FALSE) {
    ierr = grid.rescale(grid.p->Lx, grid.p->Ly, 1000.0); CHKERRQ(ierr);
  }

  // allocate observed velocity space
  ierr = createROSSVecs(); CHKERRQ(ierr);

  // update surface elev: APPLY FLOATING CRITERION EVERYWHERE to get smooth surface;
  // does this need to be done *after* afterInitHook()?
  ierr = VecCopy(vH,vh); CHKERRQ(ierr);
  ierr = VecScale(vh, 1.0 - ice.rho / ocean.rho ); CHKERRQ(ierr);
  
  // fill in temperatures at depth according to special rule:
  // temp in column equals temp at surface
  ierr = fillinTemps();  CHKERRQ(ierr);

  // zeros out vuvbar; not SIA velocities will not be computed so this will stay
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);

  ierr = verbPrintf(5,grid.com,"  [using Schoof regularization constant = %10.5e]\n",
              PetscSqr(regularizingVelocitySchoof/regularizingLengthSchoof)); CHKERRQ(ierr);

  ierr = afterInitHook(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running EISMINT ROSS ...\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::diagnosticRun() {
  PetscErrorCode  ierr;
  PetscReal       tuneparam[20];
  PetscTruth      tuneSet;
  PetscInt        numtuneparam=3;
  
  // user might enter
  //     pisms -ross -tune 1.7e8,0.1e8,2.5e8
  // [note no spaces between parameters!]
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-tune", tuneparam, &numtuneparam,
                                  &tuneSet); CHKERRQ(ierr);
  if ((tuneSet == PETSC_TRUE) && (numtuneparam != 3)) {
    SETERRQ(1,"option -tune not accompanied by three parameters");
  }

  if (tuneSet == PETSC_FALSE) {
    ierr = IceModel::diagnosticRun(); CHKERRQ(ierr);
  } else {
    // in Matlab notation:  hard = tuneparam[0]:tuneparam[1]:tuneparam[2]
    for (PetscScalar hard = tuneparam[0]; hard < tuneparam[2]+0.1*tuneparam[1];
                     hard += tuneparam[1]) { 
      ierr = verbPrintf(2,grid.com,"computing velocities with hardness (=bar B) = %10.5e ...\n",
                        hard); CHKERRQ(ierr);
      constantHardnessForMacAyeal = hard;
      ierr = IceModel::diagnosticRun(); CHKERRQ(ierr);
      ierr = finishROSS(); CHKERRQ(ierr);
      ierr = createROSSVecs(); CHKERRQ(ierr);
    }
  }
  return 0;
}


PetscErrorCode IceROSSModel::finishROSS() {
  PetscErrorCode  ierr;

  PetscTruth  ssaBCset;
  char        ssaBCfile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-ssaBC", ssaBCfile,
                               PETSC_MAX_PATH_LEN, &ssaBCset); CHKERRQ(ierr);
  if (ssaBCset == PETSC_FALSE)
    return 0;
    
  ierr = verbPrintf(2,grid.com, "\nreading EISMINT ROSS observed velocities from %s...\n",
                    ssaBCfile); CHKERRQ(ierr);
  ierr = readObservedVels(ssaBCfile); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "computing error relative to EISMINT ROSS observed velocities ...\n");
      CHKERRQ(ierr);
  ierr = computeErrorsInAccurateRegion(); CHKERRQ(ierr);
  //  ierr = readRIGGSandCompare(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "showing EISMINT ROSS observed velocities"); CHKERRQ(ierr);
  ierr = putObservedVelsCartesian(); CHKERRQ(ierr);
  ierr = updateViewers(); CHKERRQ(ierr);
  Vec vNu[2] = {vWork2d[0], vWork2d[1]};
  ierr = computeEffectiveViscosity(vNu, macayealEpsilon); CHKERRQ(ierr);
  ierr = updateNuViewers(vNu,vNu,false); CHKERRQ(ierr);
  
  PetscInt    pause_time = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, PETSC_NULL); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"\npausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }
  
  ierr = destroyROSSVecs(); CHKERRQ(ierr);
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


PetscErrorCode IceROSSModel::readObservedVels(const char *fname) {
  PetscErrorCode  ierr;

  // determine if variables exist in file and get missing values; broadcast
  int    ncid, stat, magid, aziid, accid, magExists=0, aziExists=0, accExists=0;
  double magMiss, aziMiss; 
  int    accMiss;
  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "accur", &accid); accExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "mag_obs", &magid); magExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "azi_obs", &aziid); aziExists = stat == NC_NOERR;
    stat = nc_get_att_int(ncid, accid, "missing_value",&accMiss); CHKERRQ(nc_check(stat));
    stat = nc_get_att_double(ncid, magid, "missing_value",
                             &magMiss); CHKERRQ(nc_check(stat));
    stat = nc_get_att_double(ncid, aziid, "missing_value",
                             &aziMiss); CHKERRQ(nc_check(stat));
  }
  ierr = MPI_Bcast(&accExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&magExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&aziExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&accMiss, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&magMiss, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&aziMiss, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);

  // compare comment in bootstrapFromFile_netCDF 
  Vec vzero;
  VecScatter ctx;
  ierr = VecScatterCreateToZero(g2, &ctx, &vzero); CHKERRQ(ierr);  
  ierr = getIndZero(grid.da2, g2, vzero, ctx); CHKERRQ(ierr);
  // compare comment in bootstrapFromFile_netCDF 
  if (accExists) {
    MaskInterp masklevs;
    masklevs.number_allowed = 3;
    masklevs.allowed_levels[0] = accMiss;
    masklevs.allowed_levels[1] = 0;
    masklevs.allowed_levels[2] = 1;
    ierr = ncVarToDAVec(ncid, accid, grid.da2, obsAccurate, g2, vzero, masklevs); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"accur does not exist but need to call ncVarToDAVec on it");
  }
  if (magExists) {
    ierr = ncVarToDAVec(ncid, magid, grid.da2, obsMagnitude, g2, vzero); CHKERRQ(ierr);
  } else {
    SETERRQ(2,"mag does not exist but need to call ncVarToDAVec on it");
  }
  if (aziExists) {
    ierr = ncVarToDAVec(ncid, aziid, grid.da2, obsAzimuth, g2, vzero); CHKERRQ(ierr);
  } else {
    SETERRQ(3,"azi does not exist but need to call ncVarToDAVec on it");
  }
  
  // done with scatter context
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);
  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }
  return 0;
}


PetscErrorCode IceROSSModel::computeErrorsInAccurateRegion() {
  // average over grid, where observed velocity is accurate *according to
  // 111by147grid.dat*, the difference between computed and observed u,v
  PetscErrorCode  ierr;
  PetscScalar  uerr=0.0, verr=0.0, relvecerr=0.0, accN=0.0, 
               accArea=0.0, maxcComputed=0.0, vecErrAcc = 0.0;
  PetscScalar  guerr, gverr, grelvecerr, gaccN, 
               gaccArea, gmaxcComputed, gvecErrAcc;
  PetscScalar  **azi, **mag, **acc, **ubar, **vbar, **H, **mask;
  
  const PetscScalar pi = 3.14159265358979, area = grid.p->dx * grid.p->dy;
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azi); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &mag); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsAccurate, &acc); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
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
          relvecerr += (PetscSqr(Dv) + PetscSqr(Du)) / (PetscSqr(vobs) + PetscSqr(uobs));
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

  ierr = PetscGlobalMax(&maxcComputed, &gmaxcComputed, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&accN, &gaccN, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&accArea, &gaccArea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&verr, &gverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&uerr, &guerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&relvecerr, &grelvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&vecErrAcc, &gvecErrAcc, grid.com); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"maximum computed speed in ice shelf is %10.3f (m/a)\n",
             gmaxcComputed * secpera); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"ERRORS relative to observations of Ross Ice Shelf:\n");
             CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  [number of grid points in 'accurate observed area' = %d]\n",
             (int) gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  [area of 'accurate observed area' = %9.4f (km^2)]\n",
             gaccArea / 1e6); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  following are average errors computed over 'accurate observed area':\n");
             CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average error in x-comp of vel       = %9.3f (m/a)\n",
             (gverr * secpera) / gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average error in y-comp of vel       = %9.3f (m/a)\n",
             (guerr * secpera) / gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average relative error in vector vel = %9.5f\n",
             grelvecerr / gaccN); CHKERRQ(ierr);
  gvecErrAcc = secpera * sqrt(gvecErrAcc) / sqrt(gaccArea);
  ierr = verbPrintf(2,grid.com,
             "  rms average error in vector vel      = %9.3f (m/a)\n",
             gvecErrAcc); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::putObservedVelsCartesian() {
  PetscErrorCode  ierr;
  PetscScalar **azi, **mag, **acc, **ubar, **vbar;   
  const PetscScalar pi = 3.14159265358979;

  ierr = DAVecGetArray(grid.da2, obsAzimuth, &azi); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsMagnitude, &mag); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, obsAccurate, &acc); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
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
  return 0;
}


#if 0
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
#endif
