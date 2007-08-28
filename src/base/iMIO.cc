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

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"


PetscErrorCode  IceModel::VecViewDA2Matlab(Vec l, PetscViewer v, const char *varname) {
  PetscErrorCode ierr;
  
  ierr=PetscObjectSetName((PetscObject) g2, varname); CHKERRQ(ierr);
  ierr = DALocalToGlobal(grid.da2, l, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecView(g2, v); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(v,"\n%s = reshape(%s,%d,%d);\n\n",
             varname,varname,grid.p->Mx,grid.p->My); CHKERRQ(ierr);
  return 0;
}


bool IceModel::hasSuffix(const char* fname, const char *suffix) const {
  int flen = strlen(fname);
  int slen = strlen(suffix);
  if (strcmp(fname + flen - slen, suffix) == 0) {
    return true;
  } else {
    return false;
  }
}


PetscErrorCode IceModel::initFromFile(const char *fname) {
  PetscErrorCode  ierr;

  if (hasSuffix(fname, ".pb") == true) {
    SETERRQ1(1,"ERROR: .pb format no longer supported; cannot initialize from file %s", fname);
  }
  
  if (hasSuffix(fname, ".nc") == false) {
    ierr = verbPrintf(1,grid.com,
       "WARNING:  Unknown file format for %s.  Trying to read as NetCDF.\n",fname); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2,grid.com,"initializing from NetCDF format file  %s  ...\n",
                     fname); CHKERRQ(ierr);
  ierr = initFromFile_netCDF(fname); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModel::setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID) {
  PetscErrorCode ierr;

  // read options about year of start, year of end, number of run years;
  // note grid.p->year has already been set from input file
  PetscScalar usrStartYear, usrEndYear, usrRunYears;
  PetscTruth ysSet, yeSet, ySet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ys", &usrStartYear, &ysSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ye", &usrEndYear, &yeSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-y", &usrRunYears, &ySet); CHKERRQ(ierr);
  if (ysSet == PETSC_TRUE) {
    // user option overwrites data
    ierr = setStartYear(usrStartYear); CHKERRQ(ierr);
    grid.p->year = usrStartYear;
  } else if (grid_p_year_VALID == PETSC_TRUE) {
    ierr = setStartYear(grid.p->year); CHKERRQ(ierr);
  } // else do nothing; defaults are set
  if (yeSet == PETSC_TRUE) {
    if (usrEndYear < startYear) {
      SETERRQ(1,
        "ERROR: -ye value less than -ys value (or input file year or default).\n"
        "PISM cannot run backward in time");
    }
    if (ySet == PETSC_TRUE) {
      ierr = verbPrintf(1,grid.com,"WARNING: -y option ignored.  -ye used instead.\n"); CHKERRQ(ierr);
    }
    ierr = setEndYear(usrEndYear); CHKERRQ(ierr);
  } else if (ySet == PETSC_TRUE) {
    ierr = setEndYear(usrRunYears + startYear); CHKERRQ(ierr);
  } else {
    ierr = setEndYear(DEFAULT_RUN_YEARS + startYear); CHKERRQ(ierr);
  }
  
  yearsStartRunEndDetermined = PETSC_TRUE;
  return 0;
}

  
PetscErrorCode  IceModel::writeFiles(const char* basename) {
  PetscErrorCode ierr = writeFiles(basename,PETSC_FALSE); CHKERRQ(ierr);
  return 0;
}

  
PetscErrorCode  IceModel::writeFiles(const char* basename, const PetscTruth forceFullDiagnostics) {
  PetscErrorCode ierr;

  char b[PETSC_MAX_PATH_LEN];
  char fmt[PETSC_MAX_PATH_LEN] = "n";
  char ncf[PETSC_MAX_PATH_LEN]; // netCDF format
  char mf[PETSC_MAX_PATH_LEN];  // Matlab format

  if (doPDD == PETSC_TRUE) { // want to save snow accumulation map, not net accumulation
    ierr = putBackSnowAccumPDD(); CHKERRQ(ierr);
  }
  
  // Use the defaults passed from the driver if not specified on command line.
  // We should leave space for a suffix and null byte
  strncpy(b, basename, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-of", fmt, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  // to write in both NetCDF and Matlab format, for instance:
  //     '-o foo -of nm' will generate foo.nc and foo.m

  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  if (strchr(fmt, 'p') != NULL) {
    ierr = verbPrintf(1, grid.com, "WARNING: .pb format no longer supported; writing .nc"); CHKERRQ(ierr);
    strcat(b,"_pb");  // will write basename_pb.nc
    strcat(fmt,"n");
  }

  if (strchr(fmt, 'n') != NULL) {
    strcpy(ncf, b);
    strcat(ncf, ".nc");
    PetscTruth userWantsFull;
    ierr = PetscOptionsHasName(PETSC_NULL, "-full3Dout", &userWantsFull); CHKERRQ(ierr);
    if ((forceFullDiagnostics == PETSC_TRUE) || (userWantsFull == PETSC_TRUE)) {
      ierr = verbPrintf(1, grid.com, 
            "Writing model state, with full 3D velocities, to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_diagnostic_netCDF(ncf); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(1, grid.com, "Writing model state to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_netCDF(ncf); CHKERRQ(ierr);
    }
  }

  if (strchr(fmt, 'm') != NULL) {
    strcpy(mf, b);
    strcat(mf, ".m");
    ierr = verbPrintf(1, grid.com, " ... dumping certain variables to Matlab file `%s'", mf); CHKERRQ(ierr);
    ierr = dumpToFile_Matlab(mf); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::dumpToFile_Matlab(const char *fname) {
  PetscErrorCode  ierr;
  PetscViewer  viewer;

  ierr = PetscViewerASCIIOpen(grid.com, fname, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\ndisp('iceModel output:')\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
    "\nclear year x y z xx yy zz ubar vbar c h H bed mask Tkd pmMaskkd Tid Tjd Sigmakd\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"year=%10.6f;\n",grid.p->year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "x=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",grid.p->Lx,grid.p->dx,grid.p->Lx);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "y=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",grid.p->Ly,grid.p->dy,grid.p->Ly);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n\n");  CHKERRQ(ierr);

  ierr = VecViewDA2Matlab(vh, viewer, "h"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vH, viewer, "H"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vbed, viewer, "bed"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vMask, viewer, "mask"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vubar, viewer, "ubar"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vvbar, viewer, "vbar"); CHKERRQ(ierr);

  ierr = VecPointwiseMult(vWork2d[0], vubar, vubar); CHKERRQ(ierr);
  ierr = VecPointwiseMult(vWork2d[1], vvbar, vvbar); CHKERRQ(ierr);
  ierr = VecAXPY(vWork2d[0], 1.0, vWork2d[1]); CHKERRQ(ierr);  // vWork2d[0] = ubar^2+vbar^2 = c^2 now
  PetscScalar     **c;
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &c); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      c[i][j] = sqrt(c[i][j]);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &c); CHKERRQ(ierr);
  ierr = VecScale(vWork2d[0],secpera); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[0], viewer, "c"); CHKERRQ(ierr);

  PetscScalar     ***u, **ukd, ***v, **vkd, ***w, **wkd;
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &ukd); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &vkd); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &wkd); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ukd[i][j] = u[i][j][kd];
      vkd[i][j] = v[i][j][kd];
      wkd[i][j] = w[i][j][kd];
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &ukd); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &vkd); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &wkd); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[0], viewer, "ukd"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[1], viewer, "vkd"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[2], viewer, "wkd"); CHKERRQ(ierr);

  PetscScalar     ***T, **T2, **H, **pmMask;
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &T2); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &pmMask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      T2[i][j] = T[i][j][kd];
      if (T2[i][j] + ice.beta_CC_grad * H[i][j] >= DEFAULT_MIN_TEMP_FOR_SLIDING) {
        pmMask[i][j] = 1.0;
      } else {
        pmMask[i][j] = 0.0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &T2); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &pmMask); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[0], viewer, "Tkd"); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[1], viewer, "pmMaskkd"); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Thomol = Tkd - (273.15 - H*8.66e-4);\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"hand1=figure;\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"imagesc(x,y,flipud(Thomol')), axis square, colorbar\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Tcmap = get(hand1,'ColorMap'); Tcmap(64,:)=[1 1 1];\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"set(hand1,'ColorMap',Tcmap)\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
       "title('temperature at z given by -kd (white = pressure-melting)')\n\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  PetscScalar     ***Sigma, **Sigma2;
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &Sigma2); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Sigma2[i][j] = Sigma[i][j][kd];
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &Sigma2); CHKERRQ(ierr);
  ierr = VecViewDA2Matlab(vWork2d[0], viewer, "Sigmakd"); CHKERRQ(ierr);

  // make slice along y-axis of T; requires nontrivial transfer from 3D DA-based array into new
  // type of 2D DA-based array, I think
  { // explicit scoping to reduce chance of conflicts and memory leaks
    DA           daslice;
    PetscInt     N,M,n,m;
    Vec          vTid, vTjd, gslice;
    PetscScalar  **Tid, **Tjd;

    ierr = DAGetInfo(grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                     PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
    ierr = DACreate2d(grid.com, DA_XPERIODIC, DA_STENCIL_BOX,
                      grid.p->Mz, grid.p->Mx, 1, m, 1, 1,
                      PETSC_NULL, PETSC_NULL, &daslice); CHKERRQ(ierr);
    ierr = DACreateLocalVector(daslice, &vTid); CHKERRQ(ierr);
    ierr = DACreateLocalVector(daslice, &vTjd); CHKERRQ(ierr);
    ierr = DACreateGlobalVector(daslice, &gslice); CHKERRQ(ierr);

    ierr = DAVecGetArray(daslice, vTjd, &Tjd); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt k=0; k < grid.p->Mz; k++) {
        Tjd[i][k] = T[i][jd][k];
      }
    }
    ierr = DAVecRestoreArray(daslice, vTjd, &Tjd); CHKERRQ(ierr);

    ierr = DAVecGetArray(daslice, vTid, &Tid); CHKERRQ(ierr);
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k < grid.p->Mz; k++) {
        Tid[j][k] = T[id][j][k];
      }
    }
    ierr = DAVecRestoreArray(daslice, vTid, &Tid); CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"z=linspace(0,%12.3f,%d);\n",grid.p->Lz,grid.p->Mz);
        CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject) gslice,"Tid"); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vTid, INSERT_VALUES, gslice); CHKERRQ(ierr);
    ierr = VecView(gslice, viewer); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nTid = reshape(Tid,%d,%d);\n\n",
              grid.p->Mz,grid.p->My); CHKERRQ(ierr);

    ierr=PetscObjectSetName((PetscObject) gslice,"Tjd"); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vTid, INSERT_VALUES, gslice); CHKERRQ(ierr);
    ierr = VecView(gslice, viewer); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nTjd = reshape(Tjd,%d,%d);\n\n",
              grid.p->Mz,grid.p->Mx); CHKERRQ(ierr);

    ierr = VecDestroy(gslice); CHKERRQ(ierr);
    ierr = VecDestroy(vTid); CHKERRQ(ierr);
    ierr = VecDestroy(vTjd); CHKERRQ(ierr);
    ierr = DADestroy(daslice); CHKERRQ(ierr);
  }
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  return 0;
}
