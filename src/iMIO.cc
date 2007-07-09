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
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"

/*
PARTS WHICH READ/WRITE .pb FILES ARE DEPRECATED!
*/

PetscErrorCode  IceModel::LVecLoad(DA da, Vec l, Vec g, PetscViewer v) {
  PetscErrorCode ierr;

  ierr = VecLoadIntoVector(v, g); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, l); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, l); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModel::LVecView(DA da, Vec l, Vec g, PetscViewer v) {
  PetscErrorCode ierr;
  
  ierr = DALocalToGlobal(da, l, INSERT_VALUES, g); CHKERRQ(ierr);
  ierr = VecView(g, v); CHKERRQ(ierr);
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
  PetscViewer viewer;
  PetscScalar runYears;

  if (hasSuffix(fname, ".nc") == true) {
    ierr = PetscPrintf(grid.com,
                       "initializing from NetCDF format file  %s  ...\n",
                       fname); CHKERRQ(ierr);
    ierr = initFromFile_netCDF(fname); CHKERRQ(ierr);
    return 0;
  } else if (hasSuffix(fname, ".pb") == false) {
    ierr = PetscPrintf(grid.com, "[Unknown file format."
                       "  Trying to read as PETSc binary.]\n"); CHKERRQ(ierr);
  }

  ierr = PetscPrintf(grid.com,
                     "initializing from PETSc binary format file `%s'  ...\n",
                     fname); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(grid.com, fname, FILE_MODE_READ, &viewer);
  CHKERRQ(ierr);
//  ierr = PetscBagLoad(viewer, &grid.bag); CHKERRQ(ierr);
//  ierr = PetscBagGetData(grid.bag, (void **)&(grid.p)); CHKERRQ(ierr);    
  if (relativeEndYear == PETSC_TRUE) {
    runYears = endYear - startYear;
  } else {
    runYears = endYear - grid.p->year;
  }
  ierr = setStartYear(grid.p->year); CHKERRQ(ierr);
  ierr = setRunYears(runYears); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vMask,    g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vh,       g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vH,       g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vbed,     g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vAccum,   g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vTs,      g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vGhf,     g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da2, vuplift,  g2, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da3, vT,       g3, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da3b, vTb,     g3b, viewer); CHKERRQ(ierr);
  ierr = LVecLoad(grid.da3, vtau,     g3, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);

  ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);  
  // FIXME: vHmelt should probably be part of saved state.  In this case the
  // best procedure would be to *check* if vHmelt was saved, and if so to load
  // it, otherwise to set it to zero and report that.  Similar behavior for
  // vuplift, vtau, others?  Similar behavior for many state variables in
  // loading .nc files!
  // Note: As of rev 87, vHmelt is read in the netCDF file.  Since we are
  // migrating to all netCDF I/O, vHmelt will not be a part of the Petsc saved
  // state.

  setConstantGrainSize(DEFAULT_GRAIN_SIZE);
  
  // We should not do this any more since age is part of the saved state.
  // setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);

  initialized_p = PETSC_TRUE;
  return 0;
}


PetscErrorCode  IceModel::writeFiles(const char* basename) {
  PetscErrorCode ierr;
  // Write Petsc Binary format by default
  ierr = writeFiles(basename, "n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::writeFiles(const char* basename, const char* formats) {
  PetscErrorCode ierr;
  char b[PETSC_MAX_PATH_LEN];
  char fmt[PETSC_MAX_PATH_LEN];
  char pf[PETSC_MAX_PATH_LEN];  // PETSc format
  char ncf[PETSC_MAX_PATH_LEN]; // netCDF format
  char mf[PETSC_MAX_PATH_LEN];  // Matlab format

  if (doPDD == PETSC_TRUE) { // want to save snow accumulation map, not net accumulation
    ierr = putBackSnowAccumPDD(); CHKERRQ(ierr);
  }
  
  // Use the defaults passed from the driver if not specified on command line.
  // We should leave space for a suffix and null byte
  strncpy(b, basename, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  strncpy(fmt, formats, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-of", fmt, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  // to write in both Petsc binary and Matlab format, for instance:
  //     '-o foo -of pm' will generate foo.pb and foo.m

  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  if (strchr(fmt, 'p') != NULL) {
    strcpy(pf, b);
    strcat(pf, ".pb");
    ierr = verbPrintf(2, grid.com, "Writing model state to file `%s'", pf); CHKERRQ(ierr);
    ierr = dumpToFile(pf); CHKERRQ(ierr);
  }

  if (strchr(fmt, 'n') != NULL) {
    strcpy(ncf, b);
    strcat(ncf, ".nc");
    ierr = verbPrintf(1, grid.com, "Writing model state to file `%s'", ncf); CHKERRQ(ierr);
    ierr = dumpToFile_netCDF(ncf); CHKERRQ(ierr);
  }

  if (strchr(fmt, 'm') != NULL) {
    strcpy(mf, b);
    strcat(mf, ".m");
    ierr = verbPrintf(1, grid.com, " ... dumping selected variables to Matlab file `%s'", mf); CHKERRQ(ierr);
    ierr = dumpToFile_Matlab(mf); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceModel::dumpToFile(const char *fname) {
  PetscErrorCode  ierr;
  PetscViewer viewer;

  ierr = PetscViewerBinaryOpen(grid.com, fname, FILE_MODE_WRITE,
                               &viewer); CHKERRQ(ierr);
  ierr = PetscViewerBinarySkipInfo(viewer); CHKERRQ(ierr);

//  ierr = PetscBagView(grid.bag, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vMask,    g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vh,       g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vH,       g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vbed,     g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vAccum,   g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vTs,      g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vGhf,     g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vuplift,  g2, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da3, vT,       g3, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da3b, vTb,     g3b, viewer); CHKERRQ(ierr);
  ierr = LVecView(grid.da3, vtau,     g3, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  
  // try to remove the irritating (and unused) .info file!
  char info_filename[PETSC_MAX_PATH_LEN+5];
  strcpy(info_filename, fname);
  strcat(info_filename, ".info");
  remove(info_filename);
//  if (remove(info_filename) != 0) {  // from <cstdio>; hopefully will not cause problems to use
//    SETERRQ1(1,"error deleting file %s",info_filename);
//  }

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

  ierr=PetscObjectSetName((PetscObject) g2,"h"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vh, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nh = reshape(h,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"H"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vH, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nH = reshape(H,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"bed"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vbed, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbed = reshape(bed,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"mask"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vMask, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nmask = reshape(mask,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"ubar"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vubar, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nubar = reshape(ubar,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"vbar"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vvbar, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nvbar = reshape(vbar,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) g2,"c"); CHKERRQ(ierr);  // vert-integrated hor speed in m/a
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
  ierr = LVecView(grid.da2, vWork2d[0], g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nc = reshape(c,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

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

  ierr=PetscObjectSetName((PetscObject) g2,"ukd"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vWork2d[0],  g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nukd = reshape(ukd,%d,%d);\n\n",grid.p->Mx,grid.p->My);
    CHKERRQ(ierr);
  ierr=PetscObjectSetName((PetscObject) g2,"vkd"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vWork2d[1],  g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nvkd = reshape(vkd,%d,%d);\n\n",grid.p->Mx,grid.p->My);
    CHKERRQ(ierr);
  ierr=PetscObjectSetName((PetscObject) g2,"wkd"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vWork2d[2],  g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nwkd = reshape(wkd,%d,%d);\n\n",grid.p->Mx,grid.p->My);
    CHKERRQ(ierr);


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

  ierr=PetscObjectSetName((PetscObject) g2,"Tkd"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vWork2d[0],  g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nTkd = reshape(Tkd,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"pmMaskkd"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vWork2d[1],  g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\npmMaskkd = reshape(pmMaskkd,%d,%d);\n\n",grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Thomol = Tkd - (273.15 - H*8.66e-4);\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"hand1=figure;\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"imagesc(x,y,flipud(Thomol')), axis square, colorbar\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Tcmap = get(hand1,'ColorMap'); Tcmap(64,:)=[1 1 1];\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"set(hand1,'ColorMap',Tcmap)\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
       "title('temperature at z given by -kd (white = pressure-melting)')\n\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

/*  DANGEROUS BECAUSE MATLAB USES LOTS OF MEMORY FOR THIS CONTOUR MAP
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"figure\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"contour(x,y,Tkd,200:2:300)\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"title('temperature at z given by -kd')\n\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);
*/

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
  ierr=PetscObjectSetName((PetscObject) g2,"Sigmakd"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vWork2d[0],  g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nSigmakd = reshape(Sigmakd,%d,%d);\n\n",
        grid.p->Mx,grid.p->My);
  CHKERRQ(ierr);

/*
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"figure\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"mesh(x,y,Sigmakd), colormap cool\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"title('strain heating \\Sigma at z given by -kd')\n\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);
*/

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

//    ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"figure\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"z=linspace(0,%12.3f,%d);\n",grid.p->Lz,grid.p->Mz);
        CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

    ierr=PetscObjectSetName((PetscObject) gslice,"Tid"); CHKERRQ(ierr);
    ierr = LVecView(daslice, vTid, gslice, viewer); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nTid = reshape(Tid,%d,%d);\n\n",grid.p->Mz,grid.p->My);
        CHKERRQ(ierr);

//    ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"subplot(2,1,1), surf(y,z,Tid,'LineStyle','none')\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"title('slice of temp at x given by -id')\n\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"view(2), colorbar\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

    ierr=PetscObjectSetName((PetscObject) gslice,"Tjd"); CHKERRQ(ierr);
    ierr = LVecView(daslice, vTjd, gslice, viewer); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\nTjd = reshape(Tjd,%d,%d);\n\n",grid.p->Mz,grid.p->Mx);
        CHKERRQ(ierr);

//    ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"subplot(2,1,2), surf(x,z,Tjd,'LineStyle','none')\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"title('slice of temp at y given by -jd')\n\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"view(2), colorbar\n");  CHKERRQ(ierr);
//    ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

    ierr = VecDestroy(gslice); CHKERRQ(ierr);
    ierr = VecDestroy(vTid); CHKERRQ(ierr);
    ierr = VecDestroy(vTjd); CHKERRQ(ierr);
    ierr = DADestroy(daslice); CHKERRQ(ierr);
  }
  
  // finally restore T
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

/*
  // next block shows full 3D view of T; dangerously large in practice
  ierr = PetscViewerASCIIPrintf(viewer,"\nclear z zz T\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"z=0.0:%12.3f:%12.3f;\n",grid.p->dz,grid.p->Lz);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy,zz]=meshgrid(x,y,z);\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) g3,"T"); CHKERRQ(ierr);
  ierr = LVecView(grid.da3, vT,       g3, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo on");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nT = reshape(T,%d,%d,%d);\n",
                                grid.p->Mz,grid.p->Mx,grid.p->My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"T=permute(T,[2 3 1]);\n");  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"figure(2), clf\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "sliceh=slice(xx,yy,zz,T,[x(%d)],[],[z(1) z(%d) z(%d)]);\n",
                                (grid.p->Mx+1)/2,(int) floor(0.25*static_cast<PetscScalar>(grid.p->Mz)),
                                (int) floor(0.5*static_cast<PetscScalar>(grid.p->Mz)));
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"set(sliceh,'FaceColor','interp','EdgeColor','none')\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "alpha(1.0), colorbar, title('temperature slices; surface elevation mesh'), hold on\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"xlabel('x (km)'), ylabel('y (km)'), zlabel('z (m)');\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "surfh=mesh(x,y,h); set(surfh,'EdgeColor','k','FaceAlpha',0.0), hold off\n");
  CHKERRQ(ierr);
*/

  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  return 0;
}
