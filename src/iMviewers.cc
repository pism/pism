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

PetscErrorCode IceModel::createOneViewerIfDesired(PetscViewer *viewer, 
                                   char name, const char* title) {
  PetscErrorCode ierr;
  const int SIZE = 320, bigSIZE = 600;
  PetscScalar  yTOx = (PetscScalar)grid.p->My / (PetscScalar)grid.p->Mx;
  int size, x_dim, y_dim;
  if (strchr(diagnosticBIG, name) != NULL) {
    size = bigSIZE;
  } else if (strchr(diagnostic, name) != NULL) {
    size = SIZE;
  } else {
    *viewer = PETSC_NULL;
    return 0;
  }
  if (grid.p->My > grid.p->Mx) {
    x_dim = size; y_dim = (PetscInt) ((PetscScalar) size * yTOx); 
  } else {
    y_dim = size; x_dim = (PetscInt) ((PetscScalar) size / yTOx);
  } 
  // note we reverse x_dim <-> y_dim; see IceGrid::createDA() for original reversal
  ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, title,
           PETSC_DECIDE, PETSC_DECIDE, y_dim, x_dim, viewer);  CHKERRQ(ierr);
  
  // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
  // there is a no-titles bug
  PetscDraw draw;
  ierr = PetscViewerDrawGetDraw(*viewer,0,&draw); CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw,title); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::createViewers() {
  // It is important that the createVecs() has been called before we call this.
  PetscErrorCode ierr;

  ierr = createOneViewerIfDesired(&surfHorSpeedView, '0',"hor. speed at surface (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&surfuView, '1',"u at surface (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&surfvView, '2',"v at surface (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&surfwView, '3',"w at surface (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&accumView, 'a',"M (surface accum rate; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&bedView, 'b',"b (bed elev; m above sea level)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&betaView, 'B',"log(beta) (drag coeff; log_10(Pa s m^-1))");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&speedView, 'c',"log(speed) (log_10(m/a))");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&taucView, 'C',"tau_c (till yield stress; bar=10^5Pa)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&diffusView, 'D',"D (diffusivity; m^2/s)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&tauView, 'e',"age of ice (years) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&tauMapView, 'E',"age of ice (years) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&dhView, 'f',"thickening rate dH/dt (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&ghfView, 'F',"geothermal heat flux (mW/m^2)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&gsView, 'g',"grain size (mm) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&gsMapView, 'G',"grain size (mm) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&HView, 'H',"H (thickness; m)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&hView, 'h',"h (surface elev; m above sea level)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&nuView[0], 'i',"nu*H (I offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&nuView[1], 'j',"nu*H (J offset)");  CHKERRQ(ierr);
 
  if (strchr(diagnostic, 'k') != NULL) {
    //ierr = KSPLGMonitorCreate(PETSC_NULL, "KSP Monitor", PETSC_DECIDE, PETSC_DECIDE,
    //                          PETSC_DECIDE, PETSC_DECIDE, &kspLG); CHKERRQ(ierr);
    ierr = KSPMonitorLGCreate(PETSC_NULL, "KSP Monitor", PETSC_DECIDE, PETSC_DECIDE,
                              PETSC_DECIDE, PETSC_DECIDE, &kspLG); CHKERRQ(ierr);
    //ierr = KSPSetMonitor(MacayealKSP, KSPLGMonitor, kspLG, 0); CHKERRQ(ierr);
    ierr = KSPMonitorSet(MacayealKSP, KSPMonitorLG, kspLG, 0); CHKERRQ(ierr);
  } else kspLG = PETSC_NULL;

  ierr = createOneViewerIfDesired(&basalmeltView, 'l',"basal melt rate (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&HmeltView, 'L',"basal melt water thickness (m)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&maskView, 'm',"mask (1=SHEET, 2=DRAG, 3=FLOAT)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&lognuView, 'n',"log_10(nu*H)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&NuView[0], 'N',"(nu*H)_t (I offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&NuView[1], 'N',"(nu*H)_t (J offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&upliftView, 'p',"bed uplift rate (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&slidespeedView, 'q',"log(basal sliding speed) (log_10(m/a))");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&TsView, 'r',"suRface temperature (K)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&RbView, 'R',"basal frictional heating (mW/m^2)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&SigmaView, 's',"Sigma (strain heating; K/a) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&SigmaMapView, 'S',"Sigma (strain heating; K/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&TView, 't',"T (temperature; K) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&T2View, 'T',"T (temperature; K) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&uvbarView[0], 'U',"uvbar[0] (velocity on stag grid; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&ubarView, 'u',"ubar (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&uvbarView[1], 'V',"uvbar[1] (velocity on stag grid; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vbarView, 'v',"vbar (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&uView, 'x',"u (velocity; m/a) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&umapView, 'X',"u (velocity; m/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vView, 'y',"v (velocity; m/a) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vmapView, 'Y',"v (velocity; m/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&wView, 'z',"w (velocity; m/a) at id,jd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&wmapView, 'Z',"w (velocity; m/a) at kd");  CHKERRQ(ierr);

  createViewers_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceModel::destroyViewers() {
  PetscErrorCode ierr;

  if (kspLG != PETSC_NULL) { ierr = KSPMonitorLGDestroy(kspLG); CHKERRQ(ierr); }
  if (uvbarView[0] != PETSC_NULL) { ierr = PetscViewerDestroy(uvbarView[0]); CHKERRQ(ierr); }
  if (uvbarView[1] != PETSC_NULL) { ierr = PetscViewerDestroy(uvbarView[1]); CHKERRQ(ierr); }
  if (nuView[0] != PETSC_NULL) { ierr = PetscViewerDestroy(nuView[0]); CHKERRQ(ierr); }
  if (nuView[1] != PETSC_NULL) { ierr = PetscViewerDestroy(nuView[1]); CHKERRQ(ierr); }
  if (lognuView != PETSC_NULL) { ierr = PetscViewerDestroy(lognuView); CHKERRQ(ierr); }
  if (NuView[0] != PETSC_NULL) { ierr = PetscViewerDestroy(NuView[0]); CHKERRQ(ierr); }
  if (NuView[1] != PETSC_NULL) { ierr = PetscViewerDestroy(NuView[1]); CHKERRQ(ierr); }
  if (ubarView != PETSC_NULL) { ierr = PetscViewerDestroy(ubarView); CHKERRQ(ierr); }
  if (vbarView != PETSC_NULL) { ierr = PetscViewerDestroy(vbarView); CHKERRQ(ierr); }
  if (accumView != PETSC_NULL) { ierr = PetscViewerDestroy(accumView); CHKERRQ(ierr); }
  if (bedView != PETSC_NULL) { ierr = PetscViewerDestroy(bedView); CHKERRQ(ierr); }
  if (betaView != PETSC_NULL) { ierr = PetscViewerDestroy(betaView); CHKERRQ(ierr); }
  if (taucView != PETSC_NULL) { ierr = PetscViewerDestroy(taucView); CHKERRQ(ierr); }
  if (HmeltView != PETSC_NULL) { ierr = PetscViewerDestroy(HmeltView); CHKERRQ(ierr); }
  if (basalmeltView != PETSC_NULL) { ierr = PetscViewerDestroy(basalmeltView); CHKERRQ(ierr); }
  if (ghfView != PETSC_NULL) { ierr = PetscViewerDestroy(ghfView); CHKERRQ(ierr); }
  if (upliftView != PETSC_NULL) { ierr = PetscViewerDestroy(upliftView); CHKERRQ(ierr); }
  if (maskView != PETSC_NULL) { ierr = PetscViewerDestroy(maskView); CHKERRQ(ierr); }
  if (speedView != PETSC_NULL) { ierr = PetscViewerDestroy(speedView); CHKERRQ(ierr); }
  if (slidespeedView != PETSC_NULL) { ierr = PetscViewerDestroy(slidespeedView); CHKERRQ(ierr); }
  if (HView != PETSC_NULL) { ierr = PetscViewerDestroy(HView); CHKERRQ(ierr); }
  if (hView != PETSC_NULL) { ierr = PetscViewerDestroy(hView); CHKERRQ(ierr); }
  if (diffusView != PETSC_NULL) { ierr = PetscViewerDestroy(diffusView); CHKERRQ(ierr); }
  if (dhView != PETSC_NULL) { ierr = PetscViewerDestroy(dhView); CHKERRQ(ierr); }
  if (TView != PETSC_NULL) { ierr = PetscViewerDestroy(TView); CHKERRQ(ierr); }
  if (T2View != PETSC_NULL) { ierr = PetscViewerDestroy(T2View); CHKERRQ(ierr); }
  if (TsView != PETSC_NULL) { ierr = PetscViewerDestroy(TsView); CHKERRQ(ierr); }
  if (RbView != PETSC_NULL) { ierr = PetscViewerDestroy(RbView); CHKERRQ(ierr); }
  if (uView != PETSC_NULL) { ierr = PetscViewerDestroy(uView); CHKERRQ(ierr); }
  if (vView != PETSC_NULL) { ierr = PetscViewerDestroy(vView); CHKERRQ(ierr); }
  if (wView != PETSC_NULL) { ierr = PetscViewerDestroy(wView); CHKERRQ(ierr); }
  if (surfHorSpeedView != PETSC_NULL) { ierr = PetscViewerDestroy(surfHorSpeedView); CHKERRQ(ierr); }
  if (surfuView != PETSC_NULL) { ierr = PetscViewerDestroy(surfuView); CHKERRQ(ierr); }
  if (surfvView != PETSC_NULL) { ierr = PetscViewerDestroy(surfvView); CHKERRQ(ierr); }
  if (surfwView != PETSC_NULL) { ierr = PetscViewerDestroy(surfwView); CHKERRQ(ierr); }
  if (umapView != PETSC_NULL) { ierr = PetscViewerDestroy(umapView); CHKERRQ(ierr); }
  if (vmapView != PETSC_NULL) { ierr = PetscViewerDestroy(vmapView); CHKERRQ(ierr); }
  if (wmapView != PETSC_NULL) { ierr = PetscViewerDestroy(wmapView); CHKERRQ(ierr); }
  if (SigmaView != PETSC_NULL) { ierr = PetscViewerDestroy(SigmaView); CHKERRQ(ierr); }
  if (SigmaMapView != PETSC_NULL) { ierr = PetscViewerDestroy(SigmaMapView); CHKERRQ(ierr); }
  if (gsView != PETSC_NULL) { ierr = PetscViewerDestroy(gsView); CHKERRQ(ierr); }
  if (gsMapView != PETSC_NULL) { ierr = PetscViewerDestroy(gsMapView); CHKERRQ(ierr); }
  if (tauView != PETSC_NULL) { ierr = PetscViewerDestroy(tauView); CHKERRQ(ierr); }
  if (tauMapView != PETSC_NULL) { ierr = PetscViewerDestroy(tauMapView); CHKERRQ(ierr); }

  return 0;
}


PetscErrorCode IceModel::setSoundingFromOptions() {
  PetscErrorCode ierr;
  PetscInt myid, myjd, mykd;
  PetscTruth idSet, jdSet, kdSet;
  
  myid = (grid.p->Mx - 1)/2;
  myjd = (grid.p->My - 1)/2;
  mykd = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-id", &myid, &idSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL, "-jd", &myjd, &jdSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL, "-kd", &mykd, &kdSet); CHKERRQ(ierr);
  id = myid;
  jd = myjd;
  kd = mykd;
  //  ierr = PetscPrintf(grid.com, "    !!! id,jd,kd = %3d,%3d,%3d !!!\n",id,jd,kd); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::initSounding() {
  // setup for diagnostic "sounding" views
  PetscErrorCode  ierr;

  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mbz + grid.p->Mz, &Td); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &ud); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &vd); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &wd); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &Sigmad); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &gsd); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com,PETSC_DECIDE, grid.p->Mz, &taud); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::updateSoundings() {
  PetscErrorCode ierr;
  PetscInt   Mzsum = grid.p->Mbz + grid.p->Mz;
  PetscInt   *row;

  // row gives indices only
  row = new PetscInt[Mzsum];
  for (PetscInt k=0; k < Mzsum; k++)   row[k] = k;
  
  // transfer data in [id][jd] column to soundings Vec
  if (id>=grid.xs && id<grid.xs+grid.xm && jd>=grid.ys && jd<grid.ys+grid.ym) {
    if (TView != PETSC_NULL) {
      PetscScalar ***T, ***Tb;
      ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
      ierr = VecSetValues(Td, grid.p->Mbz, row, &Tb[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValues(Td, grid.p->Mz, &row[grid.p->Mbz], &T[id][jd][0], INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
    }
    if (uView != PETSC_NULL) {
      PetscScalar ***u;
      ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
      ierr = VecSetValues(ud, grid.p->Mz, row, &u[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
    }
    if (vView != PETSC_NULL) {
      PetscScalar ***v;
      ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
      ierr = VecSetValues(vd, grid.p->Mz, row, &v[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
    }
    if (wView != PETSC_NULL) {
      PetscScalar ***w;
      ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
      ierr = VecSetValues(wd, grid.p->Mz, row, &w[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
    }
    if (SigmaView != PETSC_NULL) {
      PetscScalar ***Sigma;
      ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
      ierr = VecSetValues(Sigmad, grid.p->Mz, row, &Sigma[id][jd][0], INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    }
    if (gsView != PETSC_NULL) {
      PetscScalar ***gs;
      ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
      ierr = VecSetValues(gsd, grid.p->Mz, row, &gs[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
    }
    if (tauView != PETSC_NULL) {
      PetscScalar ***Tau;
      ierr = DAVecGetArray(grid.da3, vtau, &Tau); CHKERRQ(ierr);
      ierr = VecSetValues(gsd, grid.p->Mz, row, &Tau[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vtau, &Tau); CHKERRQ(ierr);
    }
  }
  delete [] row;  // done with setting up soundings ...

  // actually view soundings:  
  if (TView != PETSC_NULL) {
    ierr = VecAssemblyBegin(Td); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Td); CHKERRQ(ierr);
    ierr = VecView(Td, TView); CHKERRQ(ierr);
  }
  if (uView != PETSC_NULL) {
    ierr = VecAssemblyBegin(ud); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ud); CHKERRQ(ierr);
    ierr = VecScale(ud, secpera); CHKERRQ(ierr);
    ierr = VecView(ud, uView); CHKERRQ(ierr);
  }
  if (vView != PETSC_NULL) {
    ierr = VecAssemblyBegin(vd); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vd); CHKERRQ(ierr);
    ierr = VecScale(vd, secpera); CHKERRQ(ierr);
    ierr = VecView(vd, vView); CHKERRQ(ierr);
  }
  if (wView != PETSC_NULL) {
    ierr = VecAssemblyBegin(wd); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(wd); CHKERRQ(ierr);
    ierr = VecScale(wd, secpera); CHKERRQ(ierr);
    ierr = VecView(wd, wView); CHKERRQ(ierr);
  }
  if (SigmaView != PETSC_NULL) {
    ierr = VecAssemblyBegin(Sigmad); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Sigmad); CHKERRQ(ierr);
    ierr = VecScale(Sigmad, secpera); CHKERRQ(ierr);
    ierr = VecView(Sigmad, SigmaView); CHKERRQ(ierr);
  }
  if (gsView != PETSC_NULL) {
    ierr = VecAssemblyBegin(gsd); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(gsd); CHKERRQ(ierr);
    ierr = VecScale(gsd, 1.0e3); CHKERRQ(ierr); // Display in mm
    ierr = VecView(gsd, gsView); CHKERRQ(ierr);
  }
  if (tauView != PETSC_NULL) {
    ierr = VecAssemblyBegin(taud); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(taud); CHKERRQ(ierr);
    ierr = VecScale(taud, 1.0/secpera); CHKERRQ(ierr); // Display in mm
    ierr = VecView(taud, tauView); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::updateViewers() {
  // others updated elsewhere:
  // see IceModel::massBalExplicitStep() in iceModel.cc for  dhView  ("-d g")
  // see IceModel::computeMaxDiffusivity() in iMutil.cc for  diffusView ("-d f")
  // see IceModel::velocityMacAyeal() in iMmacayeal.cc for   nuView  ("-d i" or "-d j")
  //                                                   and   lognuView  ("-d n")
  //                                                   and   NuView  ("-d N")
  // see iceCompModel.cc for compensatory Sigma viewer (and redo of Sigma viewer) "-d PS"

  PetscErrorCode ierr;

  // start by updating soundings:
  ierr = updateSoundings(); CHKERRQ(ierr);
  if (T2View != PETSC_NULL) {
    ierr = getHorSliceOf3D(vT, g2, kd);
    ierr = VecView(g2, T2View); CHKERRQ(ierr);
  }
  if (tauMapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vtau, g2, kd);
    ierr = VecScale(g2, 1.0/secpera); CHKERRQ(ierr); // Display in mm
    ierr = VecView(g2, tauMapView); CHKERRQ(ierr);
  }
  if (SigmaMapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vSigma, g2, kd);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, SigmaMapView); CHKERRQ(ierr);
  }
  if (gsMapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vgs, g2, kd);
    ierr = VecScale(g2, 1.0e3); CHKERRQ(ierr); // Display in mm
    ierr = VecView(g2, gsMapView); CHKERRQ(ierr);
  }
  if (uvbarView[0] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vuvbar[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, uvbarView[0]); CHKERRQ(ierr);
  }
  if (uvbarView[1] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vuvbar[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, uvbarView[1]); CHKERRQ(ierr);
  }
  if (ubarView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vubar, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, ubarView); CHKERRQ(ierr);
  }
  if (vbarView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vvbar, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, vbarView); CHKERRQ(ierr);
  }
  if (speedView != PETSC_NULL) {
    ierr = VecPointwiseMult(vWork2d[0], vubar, vubar); CHKERRQ(ierr);
    ierr = VecPointwiseMult(vWork2d[1], vvbar, vvbar); CHKERRQ(ierr);
    ierr = VecAXPY(vWork2d[0], 1, vWork2d[1]); CHKERRQ(ierr);
    PetscScalar **a, **H;
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (H[i][j] > 0.0) {
          const PetscScalar cmpera = secpera * sqrt(a[i][j]);
          if (cmpera > 1.0e-6) {
            a[i][j] = log10(cmpera);
          } else {
            a[i][j] = -3.0;  // essentially stopped ice
          }
        } else {
          a[i][j] = -3.0; // no ice at location
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, speedView); CHKERRQ(ierr);
  }
  if (slidespeedView != PETSC_NULL) {
    PetscScalar **a, **H, **ub, **vb;
    ierr = DAVecGetArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (H[i][j] > 0.0) {
          const PetscScalar mpera = secpera * 
                     sqrt(PetscSqr(ub[i][j]) + PetscSqr(vb[i][j]));
//          if (mpera > 1.0e-6) {
          if (mpera > 1.0e-9) {
            a[i][j] = log10(mpera);
          } else {
            a[i][j] = -9.0;  // essentially stopped ice
          }
        } else {  // no ice at location
          a[i][j] = -3.0;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = VecView(g2, slidespeedView); CHKERRQ(ierr);
  }
  if (betaView != PETSC_NULL) {
    PetscScalar   **a, **beta;
    ierr = DAVecGetArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (beta[i][j] > 1.0e5) {
          a[i][j] = log10(beta[i][j]);
        } else {
          a[i][j] = 5.0;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = VecView(g2, betaView); CHKERRQ(ierr);
  }
  if (umapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vu, g2, kd);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, umapView); CHKERRQ(ierr);
  }
  if (vmapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vv, g2, kd);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, vmapView); CHKERRQ(ierr);
  }
  if (wmapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vw, g2, kd);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, wmapView); CHKERRQ(ierr);
  }
  if (surfHorSpeedView != PETSC_NULL) {
    ierr = getSurfaceValuesOf3D(vu,vWork2d[0]); CHKERRQ(ierr);
    ierr = getSurfaceValuesOf3D(vv,vWork2d[1]); CHKERRQ(ierr);
    PetscScalar **a, **us, **vs;
    ierr = DAVecGetArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &us); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[1], &vs); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        a[i][j] = sqrt(PetscSqr(us[i][j]) + PetscSqr(vs[i][j]));
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &us); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &vs); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, surfHorSpeedView); CHKERRQ(ierr);
  }
  if (surfuView != PETSC_NULL) {
    ierr = getSurfaceValuesOf3D(vu,g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, surfuView); CHKERRQ(ierr);
  }
  if (surfvView != PETSC_NULL) {
    ierr = getSurfaceValuesOf3D(vv,g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, surfvView); CHKERRQ(ierr);
  }
  if (surfwView != PETSC_NULL) {
    ierr = getSurfaceValuesOf3D(vw,g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, surfwView); CHKERRQ(ierr);
  }
  if (HView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vH, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, HView); CHKERRQ(ierr);
  }
  if (hView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vh, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, hView); CHKERRQ(ierr);
  }
  if (TsView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vTs, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, TsView); CHKERRQ(ierr);
  }
  if (accumView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vAccum, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // Display in m/a
    ierr = VecView(g2, accumView); CHKERRQ(ierr);
  }
  if (bedView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vbed, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, bedView); CHKERRQ(ierr);
  }
  if (ghfView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vGhf, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, 1000); CHKERRQ(ierr); // is in W/m^2; display in mW/m^2
    ierr = VecView(g2, ghfView); CHKERRQ(ierr);
  }
  if (RbView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vRb, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, 1000); CHKERRQ(ierr); // is in W/m^2; display in mW/m^2
    ierr = VecView(g2, RbView); CHKERRQ(ierr);
  }
  if (upliftView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vuplift, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // Display in m/a
    ierr = VecView(g2, upliftView); CHKERRQ(ierr);
  }
  if (HmeltView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vHmelt, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, HmeltView); CHKERRQ(ierr);
  }
  if (taucView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vtauc, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, 0.00001); CHKERRQ(ierr); // Display in bar
    ierr = VecView(g2, taucView); CHKERRQ(ierr);
  }
  if (basalmeltView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vbasalMeltRate, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // Display in m/a
    ierr = VecView(g2, basalmeltView); CHKERRQ(ierr);
  }
  if (maskView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vMask, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, maskView); CHKERRQ(ierr);
  }

  return 0;
}
