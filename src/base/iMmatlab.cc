// Copyright (C) 2007-2008 Ed Bueler
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
  if (strchr(matlabOutVars, name) != NULL) {
    return true;
  } else {
    return false;
  }
}


PetscErrorCode IceModel::VecView_g2ToMatlab(PetscViewer v, 
                                 const char *varname, const char *shorttitle) {
  PetscErrorCode ierr;
  
  // add Matlab comment before listing, using short title
  ierr = PetscViewerASCIIPrintf(v, "\n%%%% %s = %s\n", varname, shorttitle); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) g2, varname); CHKERRQ(ierr);
  ierr = VecView(g2, v); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(v,"\n%s = reshape(%s,%d,%d);\n\n",
             varname, varname, grid.My, grid.Mx); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::write2DToMatlab(PetscViewer v, const char scName, 
                               Vec l2, // a da2 Vec
                               const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(scName)) {
    ierr = DALocalToGlobal(grid.da2, l2, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2,scale); CHKERRQ(ierr);
    ierr = VecView_g2ToMatlab(v, tn[cIndex(scName)].name, tn[cIndex(scName)].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeSliceToMatlab(PetscViewer v, const char scName, 
                                  IceModelVec3 imv3, const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(scName)) {
    ierr = imv3.needAccessToVals(); CHKERRQ(ierr);
    ierr = imv3.getHorSlice(g2, grid.zlevels[kd]); CHKERRQ(ierr);
    ierr = imv3.doneAccessToVals(); CHKERRQ(ierr);
    ierr = VecScale(g2,scale); CHKERRQ(ierr);
    ierr = VecView_g2ToMatlab(v, tn[cIndex(scName)].name, tn[cIndex(scName)].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeSurfaceValuesToMatlab(PetscViewer v, const char scName, 
                                  IceModelVec3 imv3, const PetscScalar scale) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(scName)) {
    ierr = imv3.needAccessToVals(); CHKERRQ(ierr);
    ierr = imv3.getSurfaceValuesVec2d(g2, vH); CHKERRQ(ierr);
    ierr = imv3.doneAccessToVals(); CHKERRQ(ierr);
    ierr = VecScale(g2,scale); CHKERRQ(ierr);
    ierr = VecView_g2ToMatlab(v, tn[cIndex(scName)].name, tn[cIndex(scName)].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeSpeed2DToMatlab(
                     PetscViewer v, const char scName, Vec lu, Vec lv, // two da2 Vecs
                     const PetscScalar scale, const PetscTruth doLog, 
                     const PetscScalar log_missing) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(scName)) {
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
    ierr = VecView_g2ToMatlab(v, tn[cIndex(scName)].name, tn[cIndex(scName)].title); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeSpeedSurfaceValuesToMatlab(
                   PetscViewer v, const char scName, IceModelVec3 imv3_u, IceModelVec3 imv3_v, 
                   const PetscScalar scale, const PetscTruth doLog,
                   const PetscScalar log_missing) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(scName)) {
    ierr = imv3_u.needAccessToVals(); CHKERRQ(ierr);
    ierr = imv3_v.needAccessToVals(); CHKERRQ(ierr); 
    ierr = imv3_u.getSurfaceValuesVec2d(vWork2d[2], vH); CHKERRQ(ierr);
    ierr = imv3_v.getSurfaceValuesVec2d(vWork2d[3], vH); CHKERRQ(ierr);  
    ierr = imv3_u.doneAccessToVals(); CHKERRQ(ierr);
    ierr = imv3_v.doneAccessToVals(); CHKERRQ(ierr);
    ierr = writeSpeed2DToMatlab(v, scName, vWork2d[2], vWork2d[3], 
              scale, doLog, log_missing); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeLog2DToMatlab(
                     PetscViewer v, const char scName, Vec l, // a da2 Vec
                     const PetscScalar scale, const PetscScalar thresh,
                     const PetscScalar log_missing) {
  PetscErrorCode ierr;
  
  if (matlabOutWanted(scName)) {
    PetscScalar **a, **b;
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, l, &b); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (a[i][j] > thresh) {
          b[i][j] = log10(scale * a[i][j]);
        } else {
          b[i][j] = log_missing;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &a); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, l, &b); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView_g2ToMatlab(v, tn[cIndex(scName)].name, tn[cIndex(scName)].title); 
             CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::writeSoundingToMatlab(
                     PetscViewer v, const char scName, IceModelVec3 imv3,
                     const PetscScalar scale, const PetscTruth doTandTb) {
  if (matlabOutWanted(scName)) {
    PetscErrorCode   ierr;
    PetscInt         rlen = grid.Mz;
    PetscInt         *row;
    Vec              m;
    PetscScalar      *ivals;

    // row gives indices only
    if (doTandTb == PETSC_TRUE)   rlen += grid.Mbz;
    row = new PetscInt[rlen];
    for (PetscInt k=0; k < rlen; k++)   row[k] = k;

    ierr = VecCreateMPI(grid.com,PETSC_DECIDE, rlen, &m); CHKERRQ(ierr);

    if ((id >= grid.xs) && (id < grid.xs+grid.xm) && (jd >= grid.ys) && (jd < grid.ys+grid.ym)) {
      if (doTandTb == PETSC_TRUE) {
        ierr = T3.needAccessToVals(); CHKERRQ(ierr);
        ierr = T3.getInternalColumn(id, jd, &ivals); CHKERRQ(ierr);
        ierr = Tb3.needAccessToVals(); CHKERRQ(ierr);
        PetscScalar *ibvals;
        ierr = Tb3.getInternalColumn(id, jd, &ibvals); CHKERRQ(ierr);
        ierr = VecSetValues(m, grid.Mbz, row, ibvals, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(m, grid.Mz, &row[grid.Mbz], ivals, INSERT_VALUES);
                 CHKERRQ(ierr);
        ierr = T3.doneAccessToVals(); CHKERRQ(ierr);
        ierr = Tb3.doneAccessToVals(); CHKERRQ(ierr);
      } else {
        ierr = imv3.needAccessToVals(); CHKERRQ(ierr);
        ierr = imv3.getInternalColumn(id, jd, &ivals); CHKERRQ(ierr);
        ierr = VecSetValues(m, rlen, row, ivals, INSERT_VALUES); CHKERRQ(ierr);
        ierr = imv3.doneAccessToVals(); CHKERRQ(ierr);
      }
    }

    delete [] row;

    ierr = VecAssemblyBegin(m); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(m); CHKERRQ(ierr);
    ierr = VecScale(m,scale); CHKERRQ(ierr);

    // add Matlab comment before listing, using short title
    ierr = PetscViewerASCIIPrintf(v, "\n%%%% %s = %s\n", tn[cIndex(scName)].name, 
              tn[cIndex(scName)].title); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) m, tn[cIndex(scName)].name); CHKERRQ(ierr);
    ierr = VecView(m, v); CHKERRQ(ierr);
    
    ierr = VecDestroy(m); CHKERRQ(ierr);
  }
  return 0;
}


//! Write out selected variables in Matlab \c .m format.
/*!
Writes out the independent variables \c year, \c x, and \c y.  Then writes out variables selected with option <tt>-matv</tt> using single character names.  See Appendix C of the User's Manual.

Writes these to <tt>foo.m</tt> if option <tt>-mato foo</tt> is given.
*/
PetscErrorCode IceModel::writeMatlabVars(const char *fname) {
  PetscErrorCode ierr;
  PetscViewer  viewer;

  ierr = PetscViewerASCIIOpen(grid.com, fname, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

  // these are tools to help actually *use* the Matlab output
  ierr = PetscViewerASCIIPrintf(viewer,"\nyear = %10.6f;\n",grid.year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "x = (-%12.3f:%12.3f:%12.3f)/1000.0;\n",
                                grid.Lx,grid.dx,grid.Lx);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "y = (-%12.3f:%12.3f:%12.3f)/1000.0;\n",
                                grid.Ly,grid.dy,grid.Ly);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"z = [\n");  CHKERRQ(ierr);
  for (PetscInt k = 0; k < grid.Mz; k++) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %14.5f\n",grid.zlevels[k]);  CHKERRQ(ierr);
  } 
  ierr = PetscViewerASCIIPrintf(viewer,"];\n\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"zb = [\n");  CHKERRQ(ierr);
  for (PetscInt k = 0; k < grid.Mbz; k++) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %14.5f\n",grid.zblevels[k]);  CHKERRQ(ierr);
  } 
  ierr = PetscViewerASCIIPrintf(viewer,"];\n\n");  CHKERRQ(ierr);
  

  // now write the variables if they are wanted
  ierr = writeSpeedSurfaceValuesToMatlab(viewer, '0', u3, v3, secpera, PETSC_FALSE, 0.0);
           CHKERRQ(ierr);
  ierr = writeSurfaceValuesToMatlab(viewer, '1', u3, secpera);  CHKERRQ(ierr);
  ierr = writeSurfaceValuesToMatlab(viewer, '2', v3, secpera);  CHKERRQ(ierr);
  ierr = writeSurfaceValuesToMatlab(viewer, '3', w3, secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, '4', vub, secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, '5', vvb, secpera);  CHKERRQ(ierr);

  ierr = write2DToMatlab(viewer,'A', (pddStuffCreated == PETSC_TRUE) ? vAccumSnow : vAccum, 
                         secpera); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'C', vtauc, 0.001);  CHKERRQ(ierr);
  ierr = writeSliceToMatlab(viewer, 'E', tau3, 1.0/secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'F', vGhf, 1000.0);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'H', vH, 1.0);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'I', vtillphi, 1.0); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'L', vHmelt, 1.0);  CHKERRQ(ierr);
  // ideally would be condition: is 'Q' desired?
  ierr = computeBasalDrivingStress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr);
  ierr = getMagnitudeOf2dVectorField(vWork2d[0],vWork2d[1],vWork2d[2]); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'Q', vWork2d[2], 0.001); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'R', vRb, 1000.0);  CHKERRQ(ierr);
  ierr = writeSliceToMatlab(viewer, 'S', Sigma3, secpera);  CHKERRQ(ierr);
  ierr = writeSliceToMatlab(viewer, 'T', T3, 1.0);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'U', vuvbar[0], secpera); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'V', vuvbar[1], secpera); CHKERRQ(ierr);
  ierr = writeSliceToMatlab(viewer, 'X', u3, secpera);  CHKERRQ(ierr);
  ierr = writeSliceToMatlab(viewer, 'Y', v3, secpera);  CHKERRQ(ierr);
  ierr = writeSliceToMatlab(viewer, 'Z', w3, secpera);  CHKERRQ(ierr);

  ierr = write2DToMatlab(viewer, 'a', vAccum, secpera); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'b', vbed, 1.0); CHKERRQ(ierr);
  ierr = writeSpeed2DToMatlab(viewer, 'c', vubar, vvbar, secpera, PETSC_TRUE, -3.0);
           CHKERRQ(ierr);
  ierr = writeSoundingToMatlab(viewer,'e',tau3,1.0/secpera, PETSC_FALSE);
           CHKERRQ(ierr); // Display in years
  ierr = write2DToMatlab(viewer, 'f', vdHdt, secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'h', vh, 1.0);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'l', vbasalMeltRate, secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'm', vMask, 1.0);  CHKERRQ(ierr);
// how to do log(nu H)?:
//  ierr = createOneViewerIfDesired(&lognuView, 'n',"log_10(nu*H)");  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'p', vuplift, secpera);  CHKERRQ(ierr);
  ierr = writeSpeed2DToMatlab(viewer, 'q', vub, vvb, secpera, PETSC_TRUE, -3.0);
           CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'r', vTs, 1.0);  CHKERRQ(ierr);
  ierr = writeSoundingToMatlab(viewer,'s', Sigma3, secpera, PETSC_FALSE); CHKERRQ(ierr);
  ierr = writeSoundingToMatlab(viewer,'t', T3, 1.0, PETSC_TRUE); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'u', vubar, secpera);  CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'v', vvbar, secpera);  CHKERRQ(ierr);
  ierr = writeSoundingToMatlab(viewer,'x', u3, secpera, PETSC_FALSE); CHKERRQ(ierr);
  ierr = writeSoundingToMatlab(viewer,'y', v3, secpera, PETSC_FALSE); CHKERRQ(ierr);
  ierr = writeSoundingToMatlab(viewer,'z', w3, secpera, PETSC_FALSE); CHKERRQ(ierr);

  // now give plotting advice in Matlab comments
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
           "%% to produce a 2D color map of a map-plane view like 'H' (for instance) do:\n");
           CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> imagesc(x,y,flipud(H')), axis square, colorbar\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
           "%% to produce a 1D plot of a sounding like 'e' (for instance) do:\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> plot(z,age_sounding)\n");  CHKERRQ(ierr);
  if (matlabOutWanted('t')) { // more advice for t
    ierr = PetscViewerASCIIPrintf(viewer,
             "%% to produce a 1D plot of the sounding 't' do:\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
             "%%%%  >> plot([zb(1:end-1) z],T_sounding)\n");  CHKERRQ(ierr);
  }
  if (matlabOutWanted('T')) { // more advice for T relative to pressure-melting
    if (!matlabOutWanted('H')) { // write out 
      ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);
      ierr = write2DToMatlab(viewer, 'H', vH, 1.0);  CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,
           "%% to produce a 2D color map of homologous temp with white "
           "for at pressure-melting do:\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> Thomol = T_kd - (273.15 - H*8.66e-4);\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> hand1=figure;\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> imagesc(x,y,flipud(Thomol')), axis square, colorbar\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> Tcmap = get(hand1,'ColorMap'); Tcmap(64,:)=[1 1 1];\n");  CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
           "%%%%  >> set(hand1,'ColorMap',Tcmap)\n");  CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  return 0;
}


//! Write out the linear system of equations for the last nonlinear iteration (for SSA).
/*!
Writes out the matrix, the right-hand side, and the solution vector in a \em Matlab-readable
format, a <tt>.m</tt> file.
*/
PetscErrorCode IceModel::writeSSAsystemMatlab(Vec vNu[2]) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           file_name[PETSC_MAX_PATH_LEN], yearappend[PETSC_MAX_PATH_LEN];

  strcpy(file_name,ssaMatlabFilePrefix);
  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_%.0f.m", grid.year);
  strcat(file_name,yearappend);
  ierr = verbPrintf(2, grid.com, 
             "writing Matlab-readable file for SSA system A xsoln = rhs to file `%s' ...\n",
             file_name); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name);CHKERRQ(ierr);

  // save linear system; gives system A xsoln = rhs at last (nonlinear) iteration of SSA
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% PISM SSA linear system report.\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Writes matrix A (sparse) and vectors xsoln, rhs for the linear system\n"); 
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% (A xsoln = rhs) which was solved at the last step of the nonlinear iteration.\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Also writes the year, the coordinates x,y and their gridded versions xx,yy.\n"); 
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Also thickness (H), surface elevation (h), i-offset (staggered grid) \n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% vertically-integrated viscosity (nu_0 = nu * H), and j-offset version\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% of same thing (nu_1 = nu * H).\n\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% The thickness can be plotted by\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%%   >> imagesc(x,y,flipud(H')), axis square, colorbar \n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% for example. \n");
             CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\necho off\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAStiffnessMatrix,"A"); CHKERRQ(ierr);
  ierr = MatView(SSAStiffnessMatrix, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSARHS,"rhs"); CHKERRQ(ierr);
  ierr = VecView(SSARHS, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAX,"xsoln"); CHKERRQ(ierr);
  ierr = VecView(SSAX, viewer);CHKERRQ(ierr);

  // save coordinates (for viewing, primarily)
  ierr = PetscViewerASCIIPrintf(viewer,"year=%10.6f;\n",grid.year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "x=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",grid.Lx,grid.dx,grid.Lx);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "y=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",grid.Ly,grid.dy,grid.Ly);
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n\n");  CHKERRQ(ierr);

  // also save thickness and effective viscosity
  ierr = write2DToMatlab(viewer, 'H', vH, 1.0); CHKERRQ(ierr);
  ierr = write2DToMatlab(viewer, 'h', vh, 1.0); CHKERRQ(ierr);
  ierr = DALocalToGlobal(grid.da2, vNu[0], INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecView_g2ToMatlab(viewer, "nu_0", 
            "effective viscosity times thickness (i offset)");
            CHKERRQ(ierr);
  ierr = DALocalToGlobal(grid.da2, vNu[1], INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecView_g2ToMatlab(viewer, "nu_1", 
            "effective viscosity times thickness (j offset)");
            CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  return 0;
}

