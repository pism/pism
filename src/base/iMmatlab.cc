// Copyright (C) 2007 Ed Bueler
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
#include "iceModel.hh"


bool IceModel::matlabOutWanted(const char name) {
  if (strchr(matlabOut, name) != NULL) {
    return true;
  } else {
    return false;
  }
}


PetscErrorCode  IceModel::VecView_g2ToMatlab(PetscViewer v, 
                                 const char *varname, const char *shorttitle) {
  PetscErrorCode ierr;
  
  // add comment before listing, using hort title
  ierr = PetscViewerASCIIPrintf(v, "\n\% %s = %s\n", varname, shorttitle); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) g2, varname); CHKERRQ(ierr);
  ierr = VecView(g2, v); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(v,"\n%s = reshape(%s,%d,%d);\n\n",
             varname, varname, grid.p->Mx, grid.p->My); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                               Vec l2, // a da2 Vec
                               const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(name)) {
    ierr = DALocalToGlobal(grid.da2, l2, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2,scale); CHKERRQ(ierr);
    const int index = int(singleCharName) - int('0');
    ierr = VecView_g2ToMatlab(v, tn[index].name, tn[index].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode writeSliceToMatlab(PetscViewer v, const char singleCharName, 
                                  Vec l3, // a da3 Vec
                                  const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(name)) {
    ierr = getHorSliceOf3D(l3, g2, kd);
    ierr = VecScale(g2,scale); CHKERRQ(ierr);
    const int index = int(singleCharName) - int('0');
    ierr = VecView_g2ToMatlab(v, tn[index].name, tn[index].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode writeSurfaceValuesToMatlab(PetscViewer v, const char singleCharName, 
                                  Vec l3, // a da3 Vec
                                  const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(name)) {
    ierr = getSurfaceValuesOf3D(l3,g2); CHKERRQ(ierr);
    ierr = VecScale(g2,scale); CHKERRQ(ierr);
    const int index = int(singleCharName) - int('0');
    ierr = VecView_g2ToMatlab(v, tn[index].name, tn[index].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeMatlab(const char *fname) {
  PetscErrorCode ierr;
  PetscViewer  viewer;

  ierr = PetscViewerASCIIOpen(grid.com, fname, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

// how to do speed?
//  ierr = writeSurfaceValuesToMatlab(viewer, '0',"hor. speed at surface (m/a)");  CHKERRQ(ierr);
  ierr = writeSurfaceValuesToMatlab(viewer, '1', vu, secpera);  CHKERRQ(ierr);
  ierr = writeSurfaceValuesToMatlab(viewer, '2', vv, secpera);  CHKERRQ(ierr);
  ierr = writeSurfaceValuesToMatlab(viewer, '3', vw, secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'a', vAccum, secpera); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'b', vbed, 1.0); CHKERRQ(ierr);
// how to do speed?
//  ierr = createOneViewerIfDesired(&speedView, 'c',"log(speed) (log_10(m/a))");  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'C', );  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&tauMapView, 'E',"age of ice (years) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&dhView, 'f',"thickening rate dH/dt (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&ghfView, 'F',"geothermal heat flux (mW/m^2)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&gsMapView, 'G',"grain size (mm) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&HView, 'H',"H (thickness; m)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&hView, 'h',"h (surface elev; m above sea level)");  CHKERRQ(ierr);
 
  ierr = createOneViewerIfDesired(&basalmeltView, 'l',"basal melt rate (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&HmeltView, 'L',"basal melt water thickness (m)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&maskView, 'm',"mask (1=SHEET, 2=DRAG, 3=FLOAT)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&lognuView, 'n',"log_10(nu*H)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&upliftView, 'p',"bed uplift rate (m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&slidespeedView, 'q',"log(basal sliding speed) (log_10(m/a))");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&TsView, 'r',"suRface temperature (K)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&RbView, 'R',"basal frictional heating (mW/m^2)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&SigmaMapView, 'S',"Sigma (strain heating; K/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&T2View, 'T',"T (temperature; K) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&ubarView, 'u',"ubar (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vbarView, 'v',"vbar (velocity; m/a)");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&umapView, 'X',"u (velocity; m/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&vmapView, 'Y',"v (velocity; m/a) at kd");  CHKERRQ(ierr);
  ierr = createOneViewerIfDesired(&wmapView, 'Z',"w (velocity; m/a) at kd");  CHKERRQ(ierr);


  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  return 0;
}







PetscErrorCode IceModel::updateViewers() {
  // others updated elsewhere:
  // see IceModel::massBalExplicitStep() in iceModel.cc for  dHdtView  ("-d f")
  // see IceModel::computeMaxDiffusivity() in iMutil.cc for  diffusView ("-d D")
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
            a[i][j] = -6.0;  // essentially stopped ice
          }
        } else {
          a[i][j] = -6.0; // no ice at location
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
          if (mpera > 1.0e-6) {
            a[i][j] = log10(mpera);
          } else {
            a[i][j] = -6.0;  // essentially stopped ice
          }
        } else {  // no ice at location
          a[i][j] = -6.0;
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
    PetscScalar **a, **us, **vs, **H;
    ierr = DAVecGetArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &us); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[1], &vs); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (H[i][j] == 0.0) {
          a[i][j] = 0.0; 
        } else {
          a[i][j] = sqrt(PetscSqr(us[i][j]) + PetscSqr(vs[i][j]));
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, g2, &a); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
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
  if (dhView != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vdHdt, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr); // to report in m/a
    ierr = VecView(g2, dhView); CHKERRQ(ierr);
  }

  return 0;
}
