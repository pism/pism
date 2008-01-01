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
#include <cmath>
#include <petscda.h>
#include <petscksp.h>
#include "iceModel.hh"


PetscErrorCode IceModel::createOneViewerIfDesired(const char singleCharName) {
  return createOneViewerIfDesired(singleCharName,tn[cIndex(singleCharName)].title);
}


PetscErrorCode IceModel::createOneViewerIfDesired(const char singleCharName,
                                                  const char* title) {
  return createOneViewerIfDesired(&(runtimeViewers[cIndex(singleCharName)]),
                                  singleCharName,tn[cIndex(singleCharName)].title);
}


PetscErrorCode IceModel::getViewerDims(const PetscInt target_size, const PetscScalar Lx, const PetscScalar Ly,
                                       PetscInt *xdim, PetscInt *ydim)                             {

  // aim for smaller dimension equal to target, larger dimension larger by Ly/Lx or Lx/Ly proportion
  const double  yTOx = Ly / Lx;
  if (Ly > Lx) {
    *xdim = target_size; 
    *ydim = (PetscInt) ((double)target_size * yTOx); 
  } else {
    *ydim = target_size; 
    *xdim = (PetscInt) ((double)target_size / yTOx);
  }
  
  // if either dimension is larger than twice the target, shrink appropriately
  if (*xdim > 2 * target_size) {
    *ydim = (PetscInt) ( (double)(*ydim) * (2.0 * (double)target_size / (double)(*xdim)) );
    *xdim = 2 * target_size;
  } else if (*ydim > 2 * target_size) {
    *xdim = (PetscInt) ( (double)(*xdim) * (2.0 * (double)target_size / (double)(*ydim)) );
    *ydim = 2 * target_size;
  }
  
  // make sure minimum dimension is sufficient to see
  if (*xdim < 20)   *xdim = 20;
  if (*ydim < 20)   *ydim = 20;
  return 0;
}


PetscErrorCode IceModel::createOneViewerIfDesired(PetscViewer* v, const char singleCharName,
                                                  const char* title) {
  PetscErrorCode ierr;

  const PetscInt SIZE = 320, bigSIZE = 600;
  PetscInt size;
  if (strchr(diagnosticBIG, singleCharName) != NULL) {
    size = bigSIZE;
  } else if (strchr(diagnostic, singleCharName) != NULL) {
    size = SIZE;
  } else {
    *v = PETSC_NULL;
    return 0;
  }
  
  // viewer dims need to be determined esp. in nonsquare cases
  PetscInt x_dim, y_dim;
  ierr = getViewerDims(size, grid.p->Lx, grid.p->Ly, &x_dim, &y_dim); CHKERRQ(ierr);

  // note we reverse x_dim <-> y_dim; see IceGrid::createDA() for original reversal
  ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, title,
           PETSC_DECIDE, PETSC_DECIDE, y_dim, x_dim, v);  CHKERRQ(ierr);
  
  // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
  // there is a no titles bug
  PetscDraw draw;
  ierr = PetscViewerDrawGetDraw(*v,0,&draw); CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw,title); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::createViewers() {
  // It is important that the createVecs() has been called before we call this.
  // See iMnames.cc and IceModel::tn for most short titles on viewers.

  if (createViewers_done == PETSC_TRUE)
    return 0;

  PetscErrorCode ierr;
  const int nv = 46; // number of viewers in use
  char viewsInUse[nv] = {'0','1','2','3','4','5',
                         'A','B','C','D','E','F','G','H','L','Q','R','S',
                                 'T','U','V','X','Y','Z',
                         'a','b','c','e','f','g','h','i','j','l',
                                 'm','n','p','q','r','s','t','u','v','x',
                                 'y','z'};

  for (PetscInt nn = 0; nn < nv; nn++) {
    ierr = createOneViewerIfDesired(viewsInUse[nn]); CHKERRQ(ierr);
  } 

  if (strchr(diagnostic, 'k') != NULL) {
    ierr = KSPMonitorLGCreate(PETSC_NULL, "KSP Monitor", PETSC_DECIDE, PETSC_DECIDE,
                              PETSC_DECIDE, PETSC_DECIDE, &kspLG); CHKERRQ(ierr);
    ierr = KSPMonitorSet(SSAKSP, KSPMonitorLG, kspLG, 0); CHKERRQ(ierr);
  } else kspLG = PETSC_NULL;

//  ierr = createOneViewerIfDesired(&NuView[0], 'N',"(nu*H)_t (I offset)");  CHKERRQ(ierr);
//  ierr = createOneViewerIfDesired(&NuView[1], 'N',"(nu*H)_t (J offset)");  CHKERRQ(ierr);

  // allocate space for soundings
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mbz + grid.p->Mz - 1, &Td); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &ud); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &vd); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &wd); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &Sigmad); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &gsd); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &taud); CHKERRQ(ierr);

  createViewers_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceModel::destroyViewers() {
  PetscErrorCode ierr;

  for (PetscInt nn = 0; nn < tnN; nn++) {
    if (runtimeViewers[nn] != PETSC_NULL) {
      ierr = PetscViewerDestroy(runtimeViewers[nn]); CHKERRQ(ierr);
    }
  } 
  
  if (kspLG != PETSC_NULL) { ierr = KSPMonitorLGDestroy(kspLG); CHKERRQ(ierr); }

  if (Td != PETSC_NULL) { ierr = VecDestroy(Td); CHKERRQ(ierr); }
  if (ud != PETSC_NULL) { ierr = VecDestroy(ud); CHKERRQ(ierr); }
  if (vd != PETSC_NULL) { ierr = VecDestroy(vd); CHKERRQ(ierr); }
  if (wd != PETSC_NULL) { ierr = VecDestroy(wd); CHKERRQ(ierr); }
  if (Sigmad != PETSC_NULL) { ierr = VecDestroy(Sigmad); CHKERRQ(ierr); }
  if (gsd != PETSC_NULL) { ierr = VecDestroy(gsd); CHKERRQ(ierr); }
  if (taud != PETSC_NULL) { ierr = VecDestroy(taud); CHKERRQ(ierr); }

  return 0;
}


PetscErrorCode IceModel::updateOneSounding(
                  const char scName, Vec l, const PetscScalar scale) {
  PetscErrorCode ierr;
  if (runtimeViewers[cIndex(scName)] != PETSC_NULL) {
    ierr = VecAssemblyBegin(l); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(l); CHKERRQ(ierr);
    ierr = VecScale(l,scale); CHKERRQ(ierr);
    ierr = VecView(l, runtimeViewers[cIndex(scName)]); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::updateSoundings() {
  PetscErrorCode ierr;
  PetscInt   Mzsum = grid.p->Mbz + grid.p->Mz - 1;
  PetscInt   *row;

  // row gives indices only
  row = new PetscInt[Mzsum];
  for (PetscInt k=0; k < Mzsum; k++)   row[k] = k;
  
  // transfer data in [id][jd] column to soundings Vec
  if (id>=grid.xs && id<grid.xs+grid.xm && jd>=grid.ys && jd<grid.ys+grid.ym) {
    if (runtimeViewers[cIndex('t')] != PETSC_NULL) {
      PetscScalar ***T, ***Tb;
      ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
      // note Tb[][][Mbz-1] duplicates T[][][0] and it should not be displayed twice
      ierr = VecSetValues(Td, grid.p->Mbz - 1, row, &Tb[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValues(Td, grid.p->Mz, &row[grid.p->Mbz - 1], &T[id][jd][0], INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
    }
    if (runtimeViewers[cIndex('x')] != PETSC_NULL) {
      PetscScalar ***u;
      ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
      ierr = VecSetValues(ud, grid.p->Mz, row, &u[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
    }
    if (runtimeViewers[cIndex('y')] != PETSC_NULL) {
      PetscScalar ***v;
      ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
      ierr = VecSetValues(vd, grid.p->Mz, row, &v[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
    }
    if (runtimeViewers[cIndex('z')] != PETSC_NULL) {
      PetscScalar ***w;
      ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
      ierr = VecSetValues(wd, grid.p->Mz, row, &w[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
    }
    if (runtimeViewers[cIndex('s')] != PETSC_NULL) {
      PetscScalar ***Sigma;
      ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
      ierr = VecSetValues(Sigmad, grid.p->Mz, row, &Sigma[id][jd][0], INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    }
    if (runtimeViewers[cIndex('g')] != PETSC_NULL) {
      PetscScalar ***gs;
      ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
      ierr = VecSetValues(gsd, grid.p->Mz, row, &gs[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
    }
    if (runtimeViewers[cIndex('e')] != PETSC_NULL) {
      PetscScalar ***Tau;
      ierr = DAVecGetArray(grid.da3, vtau, &Tau); CHKERRQ(ierr);
      ierr = VecSetValues(gsd, grid.p->Mz, row, &Tau[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vtau, &Tau); CHKERRQ(ierr);
    }
  }
  delete [] row;  // done with setting up soundings ...

  // actually view soundings:  
  ierr = updateOneSounding('e',taud,1.0/secpera); CHKERRQ(ierr); // Display in years
  ierr = updateOneSounding('g',gsd,1000.0); CHKERRQ(ierr); // Display in mm
  ierr = updateOneSounding('s',Sigmad,secpera); CHKERRQ(ierr);
  ierr = updateOneSounding('t',Td,1.0); CHKERRQ(ierr);
  ierr = updateOneSounding('x',ud,secpera); CHKERRQ(ierr);
  ierr = updateOneSounding('y',vd,secpera); CHKERRQ(ierr);
  ierr = updateOneSounding('z',wd,secpera); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::update2DViewer(const char scName, Vec l2, // a da2 Vec
                                        const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (runtimeViewers[cIndex(scName)] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, l2, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, scale); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex(scName)]); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::updateSliceViewer(const char scName, Vec l3, // a da3 Vec
                                           const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (runtimeViewers[cIndex(scName)] != PETSC_NULL) {
    ierr = getHorSliceOf3D(l3, g2, kd);
    ierr = VecScale(g2, scale); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex(scName)]); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::updateSurfaceValuesViewer(const char scName, Vec l3, // a da3 Vec
                                                   const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (runtimeViewers[cIndex(scName)] != PETSC_NULL) {
    ierr = getSurfaceValuesOf3D(l3,g2); CHKERRQ(ierr);
    ierr = VecScale(g2, scale); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex(scName)]); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::updateSpeed2DViewer(
                     const char scName, Vec lu, Vec lv, // two da2 Vecs
                     const PetscScalar scale, const PetscTruth doLog, 
                     const PetscScalar log_missing) {
  PetscErrorCode ierr;
  
  if (runtimeViewers[cIndex(scName)] != PETSC_NULL) {
    PetscScalar **a, **H;

    ierr = VecPointwiseMult(vWork2d[0], lu, lu); CHKERRQ(ierr);
    ierr = VecPointwiseMult(vWork2d[1], lv, lv); CHKERRQ(ierr);
    ierr = VecAXPY(vWork2d[0], 1, vWork2d[1]); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (doLog == PETSC_TRUE) {
          if (H[i][j] > 0.0) {
            const PetscScalar cmpera = scale * sqrt(a[i][j]);
            if (cmpera > 1.0e-6) {
              a[i][j] = log10(cmpera);
            } else {
              a[i][j] = log_missing;  // essentially stopped ice
            }
          } else {
            a[i][j] = log_missing; // no ice at location
          }
        } else { // don't do log
          a[i][j] = scale * sqrt(a[i][j]);
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex(scName)]); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::updateSpeedSurfaceValuesViewer(
                   const char scName, Vec lu, Vec lv, // two da3 Vecs
                   const PetscScalar scale, const PetscTruth doLog,
                   const PetscScalar log_missing) {
  PetscErrorCode ierr;
  
  ierr = getSurfaceValuesOf3D(lu,vWork2d[2]); CHKERRQ(ierr);
  ierr = getSurfaceValuesOf3D(lv,vWork2d[3]); CHKERRQ(ierr);
  ierr = updateSpeed2DViewer(scName, vWork2d[2], vWork2d[3], 
            scale, doLog, log_missing); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::updateLog2DViewer(
                     const char scName, Vec l, // a da2 Vec
                     const PetscScalar scale, const PetscScalar thresh,
                     const PetscScalar log_missing) {
  PetscErrorCode ierr;
  
  if (runtimeViewers[cIndex(scName)] != PETSC_NULL) {
    PetscScalar **a, **b;
    ierr = DAVecGetArray(grid.da2, l, &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &b); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (a[i][j] > thresh) {
          b[i][j] = log10(scale * a[i][j]);
        } else {
          b[i][j] = log_missing;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &b); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, l, &a); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex(scName)]); CHKERRQ(ierr);
  }
  return 0;
}


//! Update the runtime graphical viewers.
/*! At every time step the graphical viewers are updated.  The user specifies these viewers
by the options <tt>-d</tt> \em list or <tt>-dbig</tt> \em list where \em list is a list of single
character names of the viewers (a list with no spaces).  See an appendix of the User's Manual
for the names.

Most viewers are updated by this routing, but some other are updated elsewhere:
  \li see computeMaxDiffusivity() in iMutil.cc for  diffusView ("-d D")
  \li see updateNuViewers() for   nuView  ("-d i" or "-d j")  and   lognuView  ("-d n")
        and   NuView  ("-d N")
  \li see iceCompModel.cc for compensatory Sigma viewer (and redo of Sigma viewer) "-d PS".
 */
PetscErrorCode IceModel::updateViewers() {

  PetscErrorCode ierr;

  // start by updating soundings:
  ierr = updateSoundings(); CHKERRQ(ierr);
  
  ierr = updateSpeedSurfaceValuesViewer('0', vu, vv, secpera, PETSC_FALSE, 0.0); CHKERRQ(ierr);
  ierr = updateSurfaceValuesViewer('1', vu, secpera); CHKERRQ(ierr);
  ierr = updateSurfaceValuesViewer('2', vv, secpera); CHKERRQ(ierr);
  ierr = updateSurfaceValuesViewer('3', vw, secpera); CHKERRQ(ierr);
  ierr = update2DViewer('4', vub, secpera); CHKERRQ(ierr);
  ierr = update2DViewer('5', vvb, secpera); CHKERRQ(ierr);

  ierr = update2DViewer('A', (pddStuffCreated == PETSC_TRUE) ? vAccumSnow : vAccum, secpera); CHKERRQ(ierr);
  ierr = updateLog2DViewer('B', vbeta, 1.0, 1.0e5, 5.0); CHKERRQ(ierr);
  ierr = update2DViewer('C', vtauc, 0.001); CHKERRQ(ierr); // display in kPa
  ierr = updateSliceViewer('E', vtau, 1.0/secpera); CHKERRQ(ierr); // display in years
  ierr = update2DViewer('F', vGhf, 1000.0); CHKERRQ(ierr); // is in W/m^2; display in mW/m^2
  ierr = updateSliceViewer('G', vgs, 1000.0); CHKERRQ(ierr); // in mm
  ierr = update2DViewer('H', vH, 1.0); CHKERRQ(ierr);
  ierr = update2DViewer('L', vHmelt, 1.0); CHKERRQ(ierr);
  ierr = computeBasalDrivingStress(vWork2d[0]); CHKERRQ(ierr);
  ierr = update2DViewer('Q', vWork2d[0], 0.001); CHKERRQ(ierr); // Display in kPa
  ierr = update2DViewer('R', vRb, 1000.0); CHKERRQ(ierr); // is in W/m^2; display in mW/m^2
  ierr = updateSliceViewer('S', vSigma, secpera); CHKERRQ(ierr);
  ierr = updateSliceViewer('T', vT, 1.0); CHKERRQ(ierr);
  ierr = update2DViewer('U', vuvbar[0], secpera); CHKERRQ(ierr);
  ierr = update2DViewer('V', vuvbar[1], secpera); CHKERRQ(ierr);
  ierr = updateSliceViewer('X', vu, secpera); CHKERRQ(ierr);
  ierr = updateSliceViewer('Y', vv, secpera); CHKERRQ(ierr);
  ierr = updateSliceViewer('Z', vw, secpera); CHKERRQ(ierr);

  ierr = update2DViewer('a', vAccum, secpera); CHKERRQ(ierr);
  ierr = update2DViewer('b', vbed, 1.0); CHKERRQ(ierr);
  ierr = updateSpeed2DViewer('c', vubar, vvbar, secpera, PETSC_TRUE, -3.0); CHKERRQ(ierr);
  // 'e' is sounding
  ierr = update2DViewer('f', vdHdt, secpera); CHKERRQ(ierr);
  // 'g' is sounding
  ierr = update2DViewer('h', vh, 1.0); CHKERRQ(ierr);
  ierr = update2DViewer('l', vbasalMeltRate, secpera); CHKERRQ(ierr);
  ierr = update2DViewer('m', vMask, 1.0); CHKERRQ(ierr);
  ierr = update2DViewer('p', vuplift, secpera); CHKERRQ(ierr);
  ierr = updateSpeed2DViewer('q', vub, vvb, secpera, PETSC_TRUE, -3.0); CHKERRQ(ierr);
  ierr = update2DViewer('r', vTs, 1.0); CHKERRQ(ierr);
  // 's' is sounding
  // 't' is sounding
  ierr = update2DViewer('u', vubar, secpera); CHKERRQ(ierr);
  ierr = update2DViewer('v', vvbar, secpera); CHKERRQ(ierr);
  // 'x' is sounding
  // 'y' is sounding
  // 'z' is sounding

  return 0;
}


PetscErrorCode IceModel::updateNuViewers(Vec vNu[2], Vec vNuOld[2], bool updateNu_tView) {
  // this one is called when solving an SSA system
  PetscErrorCode ierr;
  if (runtimeViewers[cIndex('n')] != PETSC_NULL) {
    PetscScalar  **nui, **nuj, **gg;  
    ierr = DAVecGetArray(grid.da2, g2, &gg); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vNu[0], &nui); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vNu[1], &nuj); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscReal avnu = 0.5 * (nui[i][j] + nuj[i][j]);
        if (avnu > 1.0e14) {
          gg[i][j] = log10(avnu);
        } else {
          gg[i][j] = 14.0;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vNu[0], &nui); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vNu[1], &nuj); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, g2, &gg); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex('n')]); CHKERRQ(ierr);
  }
  if (runtimeViewers[cIndex('i')] != PETSC_NULL && runtimeViewers[cIndex('j')] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vNu[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex('i')]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vNu[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex('j')]); CHKERRQ(ierr);
  }
/*
  if ((NuView[0] != PETSC_NULL) && (NuView[1] != PETSC_NULL) && updateNu_tView) {
    // note vNuOld[] contain *difference* of nu after testConvergenceofNu()
    ierr = DALocalToGlobal(grid.da2, vNuOld[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, NuView[0]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vNuOld[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, NuView[1]); CHKERRQ(ierr);
  }
*/
  return 0;
}

