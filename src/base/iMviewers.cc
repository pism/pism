// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

/*
PetscErrorCode IceModel::createViewers() {

  if ( (strchr(diagnostic, 'k') != NULL) || (strchr(diagnosticBIG, 'k') != NULL) ) {
    ierr = KSPMonitorLGCreate(PETSC_NULL, "KSP Monitor", PETSC_DECIDE, PETSC_DECIDE,
                              PETSC_DECIDE, PETSC_DECIDE, &kspLG); CHKERRQ(ierr);
    ierr = KSPMonitorSet(SSAKSP, KSPMonitorLG, kspLG, 0); CHKERRQ(ierr);
  } else kspLG = PETSC_NULL;

  return 0;
}
*/

/*
PetscErrorCode IceModel::destroyViewers() {
  PetscErrorCode ierr;
  
  if (kspLG != PETSC_NULL) { ierr = KSPMonitorLGDestroy(kspLG); CHKERRQ(ierr); }

  return 0;
}
*/

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
PetscErrorCode IceModel::update_viewers() {
  PetscErrorCode ierr;

  return 0;
}


//PetscErrorCode IceModel::updateNuViewers(IceModelVec2 vNu[2], IceModelVec2 /*vNuOld*/[2], bool /*updateNu_tView*/) {
  // this one is called when solving an SSA system
/*
  PetscErrorCode ierr;
  if (runtimeViewers[cIndex('n')] != PETSC_NULL) {
    PetscScalar  **nui, **nuj, **gg;  
    ierr = DAVecGetArray(grid.da2, g2, &gg); CHKERRQ(ierr);
    ierr = vNu[0].get_array(nui); CHKERRQ(ierr);
    ierr = vNu[1].get_array(nuj); CHKERRQ(ierr);
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
    ierr = vNu[0].end_access(); CHKERRQ(ierr);
    ierr = vNu[1].end_access(); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, g2, &gg); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex('n')]); CHKERRQ(ierr);
  }
  if (runtimeViewers[cIndex('i')] != PETSC_NULL && runtimeViewers[cIndex('j')] != PETSC_NULL) {
    ierr = vNu[0].copy_to_global(g2); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex('i')]); CHKERRQ(ierr);
    ierr = vNu[1].copy_to_global(g2); CHKERRQ(ierr);
    ierr = VecView(g2, runtimeViewers[cIndex('j')]); CHKERRQ(ierr);
  }
//   if ((NuView[0] != PETSC_NULL) && (NuView[1] != PETSC_NULL) && updateNu_tView) {
//     // note vNuOld[] contain *difference* of nu after testConvergenceofNu()
//     ierr = DALocalToGlobal(grid.da2, vNuOld[0], INSERT_VALUES, g2); CHKERRQ(ierr);
//     ierr = VecView(g2, NuView[0]); CHKERRQ(ierr);
//     ierr = DALocalToGlobal(grid.da2, vNuOld[1], INSERT_VALUES, g2); CHKERRQ(ierr);
//     ierr = VecView(g2, NuView[1]); CHKERRQ(ierr);
  }
  return 0;
}
*/
