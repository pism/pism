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
#include <cmath>
#include <petscda.h>
#include "iceModel.hh"


PetscErrorCode IceModel::createOneViewerIfDesired(PetscViewer *viewer, 
                                   char name, const char* title) {
  PetscErrorCode ierr;
  
  if (strchr(diagnosticBIG, name) != NULL) {
    const int bigsize = 600;
    ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, title,
             PETSC_DECIDE, PETSC_DECIDE, bigsize, bigsize, viewer);  CHKERRQ(ierr);
  } else if (strchr(diagnostic, name) != NULL) {
    ierr = PetscViewerDrawOpen(grid.com, PETSC_NULL, title,
             PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, viewer);  CHKERRQ(ierr);
  } else {
    *viewer = PETSC_NULL;
  }
  return 0;
}


PetscErrorCode IceModel::createViewers() {
  // It is important that the createVecs() has been called before we call this.
  PetscErrorCode ierr;

  ierr = createOneViewerIfDesired(&accumView, 'a',"M (surface accum rate; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&bedView, 'b',"b (bed elev; m above sea level)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&speedView, 'c',"log(speed) (log_10(m/a))");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&gsView, 'd',"grain size (mm)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&gsMapView, 'D',"grain size (mm) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&tauView, 'e',"age of ice (years)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&tauMapView, 'E',"age of ice (years) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&diffusView, 'f',"D (diffusivity; m^2/s)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&dhView, 'g',"thickening rate dH/dt (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&ghfView, 'G',"geothermal heat flux (mW/m^2)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&HView, 'H',"H (thickness; m)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&hView, 'h',"h (surface elev; m above sea level)");  CHKERRQ(ierr);

  if (strchr(diagnostic, 'k') != NULL) {
    ierr = KSPLGMonitorCreate(PETSC_NULL, "KSP Monitor", PETSC_DECIDE, PETSC_DECIDE,
                              PETSC_DECIDE, PETSC_DECIDE, &kspLG); CHKERRQ(ierr);
    ierr = KSPSetMonitor(MacayealKSP, KSPLGMonitor, kspLG, 0); CHKERRQ(ierr);
  } else kspLG = PETSC_NULL;

  ierr = createOneViewerIfDesired(&basalmeltView, 'l',"basal melt rate (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&HmeltView, 'L',"basal melt water thickness (m)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&maskView, 'm',"mask (1=SHEET, 2=DRAG, 3=FLOAT)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&nuView[0], 'n',"nu (I offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&nuView[1], 'n',"nu (J offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&NuView[0], 'N',"nu_t (I offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&NuView[1], 'N',"nu_t (J offset)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&upliftView, 'p',"bed uplift rate (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&slidespeedView, 'q',"log(basal sliding speed) (log_10(m/a))");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&TsView, 'r',"suRface temperature (K)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&RbView, 'R',"basal frictional heating (mW/m^2)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&SigmaView, 's',"Sigma (strain heating; K/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&SigmaMapView, 'S',"Sigma (strain heating; K/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&TView, 't',"T (temperature; K)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&T2View, 'T',"T (temperature; K) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&uvbarView[0], 'U',"uvbar[0] (velocity on stag grid; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&uvbarView[1], 'V',"uvbar[1] (velocity on stag grid; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&ubarView, 'u',"ubar (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vbarView, 'v',"vbar (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&uView, 'x',"u (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&umapView, 'X',"u (velocity; m/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vView, 'y',"v (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vmapView, 'Y',"v (velocity; m/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&wView, 'z',"w (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&wmapView, 'Z',"w (velocity; m/a) at kd");  CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::destroyViewers() {
  PetscErrorCode ierr;

  if (kspLG != PETSC_NULL) { ierr = KSPLGMonitorDestroy(kspLG); CHKERRQ(ierr); }
  if (uvbarView[0] != PETSC_NULL) { ierr = PetscViewerDestroy(uvbarView[0]); CHKERRQ(ierr); }
  if (uvbarView[1] != PETSC_NULL) { ierr = PetscViewerDestroy(uvbarView[1]); CHKERRQ(ierr); }
  if (nuView[0] != PETSC_NULL) { ierr = PetscViewerDestroy(nuView[0]); CHKERRQ(ierr); }
  if (nuView[1] != PETSC_NULL) { ierr = PetscViewerDestroy(nuView[1]); CHKERRQ(ierr); }
  if (NuView[0] != PETSC_NULL) { ierr = PetscViewerDestroy(NuView[0]); CHKERRQ(ierr); }
  if (NuView[1] != PETSC_NULL) { ierr = PetscViewerDestroy(NuView[1]); CHKERRQ(ierr); }
  if (ubarView != PETSC_NULL) { ierr = PetscViewerDestroy(ubarView); CHKERRQ(ierr); }
  if (vbarView != PETSC_NULL) { ierr = PetscViewerDestroy(vbarView); CHKERRQ(ierr); }
  if (accumView != PETSC_NULL) { ierr = PetscViewerDestroy(accumView); CHKERRQ(ierr); }
  if (bedView != PETSC_NULL) { ierr = PetscViewerDestroy(bedView); CHKERRQ(ierr); }
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


PetscErrorCode IceModel::updateViewers() {
  // see IceModel::massBalExplicitStep() in iceModel.cc for  dhView  ("-d g")
  // see IceModel::computeMaxDiffusivity() in iMutil.cc for  diffusView ("-d f")
  // see IceModel::velocityMacAyeal() in iMmacayeal.cc for   nuView  ("-d n")
  //                                                   and   NuView  ("-d N")
  // see iceCompModel.cc for compensatory Sigma viewers

  PetscErrorCode ierr;
  PetscScalar ***T, ***Tb, ***u, ***v, ***w, **ub, **vb, ***Sigma, 
              **sigma, **T2, **a, **H, ***gs, **gs2, ***Tau, **tau,
              **umap, **vmap, **wmap;
  PetscInt   Mzsum = grid.p->Mbz + grid.p->Mz;
  PetscInt   *row;

  // set up (quantity) vs z graphs in sounding:
  row = new PetscInt[Mzsum];
  for (PetscInt k=0; k < Mzsum; k++)   row[k] = k;
  if (id>=grid.xs && id<grid.xs+grid.xm && jd>=grid.ys && jd<grid.ys+grid.ym) {
    if (TView != PETSC_NULL) {
      ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
      ierr = VecSetValues(Td, grid.p->Mbz, row, &Tb[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValues(Td, grid.p->Mz, &row[grid.p->Mbz], &T[id][jd][0], INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
    }
    if (uView != PETSC_NULL) {
      ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
      ierr = VecSetValues(ud, grid.p->Mz, row, &u[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
    }
    if (vView != PETSC_NULL) {
      ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
      ierr = VecSetValues(vd, grid.p->Mz, row, &v[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
    }
    if (wView != PETSC_NULL) {
      ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
      ierr = VecSetValues(wd, grid.p->Mz, row, &w[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
    }
    if (SigmaView != PETSC_NULL) {
      ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
      ierr = VecSetValues(Sigmad, grid.p->Mz, row, &Sigma[id][jd][0], INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    }
    if (gsView != PETSC_NULL) {
      ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
      ierr = VecSetValues(gsd, grid.p->Mz, row, &gs[id][jd][0], INSERT_VALUES); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
    }
    if (tauView != PETSC_NULL) {
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
  // done with soundings

  // map-plane views:
  if (T2View != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, g2, &T2); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        T2[i][j] = T[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, g2, &T2); CHKERRQ(ierr);
    ierr = VecView(g2, T2View); CHKERRQ(ierr);
  }
  if (tauMapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &tau); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vtau, &Tau); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        tau[i][j] = Tau[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &tau); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vtau, &Tau); CHKERRQ(ierr);
    ierr = VecScale(g2, 1.0/secpera); CHKERRQ(ierr); // Display in mm
    ierr = VecView(g2, tauMapView); CHKERRQ(ierr);
  }
  if (SigmaMapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &sigma); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        sigma[i][j] = Sigma[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &sigma); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, SigmaMapView); CHKERRQ(ierr);
  }
  if (gsMapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &gs2); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        gs2[i][j] = gs[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &gs2); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
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
    ierr = DAVecGetArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (H[i][j] > 0.0) {
          const PetscScalar mpera = secpera * 
                     sqrt(PetscSqr(ub[i][j]) + PetscSqr(vb[i][j]));
          if (mpera > 1.0e-6) {
            a[i][j] = log10(mpera);
          } else {
            a[i][j] = -3.0;  // essentially stopped ice
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
  if (umapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &umap); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        umap[i][j] = u[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &umap); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, umapView); CHKERRQ(ierr);
  }
  if (vmapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &vmap); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        vmap[i][j] = v[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &vmap); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, vmapView); CHKERRQ(ierr);
  }
  if (wmapView != PETSC_NULL) {
    ierr = DAVecGetArray(grid.da2, g2, &wmap); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        wmap[i][j] = w[i][j][kd];
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &wmap); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // m/a
    ierr = VecView(g2, wmapView); CHKERRQ(ierr);
  }
  if (HView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vH, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, HView); CHKERRQ(ierr);
  }
  if (hView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vh, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, hView); CHKERRQ(ierr);
    //ierr = VecView(g2, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    //         PetscDraw hDraw;
    //         ierr = PetscViewerDrawGetDraw(hView, 0, &hDraw); CHKERRQ(ierr);
    //         ierr = PetscDrawLGSetLimits(hDraw, 0, 200, 0, 250); CHKERRQ(ierr);
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
