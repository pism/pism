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

#include <sstream>
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
  set<string>::iterator i;

  // map-plane viewers
  for (i = map_viewers.begin(); i != map_viewers.end(); ++i) {
    IceModelVec *v = variables.get(*i);

    // if not found, try to compute:
    if (v == NULL) {
      ierr = compute_by_name(*i, v); CHKERRQ(ierr);
    }

    // if still not found, ignore
    if (v == NULL)
      continue;

    GridType dims = v->grid_type();
    bool big_viewer = (big_viewers[*i]) || (big_viewers[*i + "_map"]);

    switch(dims) {
    case GRID_2D:
      {
	IceModelVec2 *v2 = dynamic_cast<IceModelVec2*>(v);
	if (v2 == NULL) SETERRQ(1,"grid_type() gives GRID_2D but dynamic_cast gives a NULL");
	ierr = v2->view(g2, big_viewer); CHKERRQ(ierr);
	break;
      }
    case GRID_3D:
      {
	IceModelVec3 *v3 = dynamic_cast<IceModelVec3*>(v);
	if (v3 == NULL) SETERRQ(1,"grid_type() gives GRID_3D but dynamic_cast gives a NULL");
	ierr = v3->view_surface(vH, g2, big_viewer); CHKERRQ(ierr);
	break;
      }
    case GRID_3D_BEDROCK:
      {
	ierr = PetscPrintf(grid.com,
			   "PISM ERROR: map-plane views of bedrock quantities are not supported.\n");
	CHKERRQ(ierr);
	PetscEnd();
      }
    }
  }

  // slice viewers:
  for (i = slice_viewers.begin(); i != slice_viewers.end(); ++i) {
    IceModelVec *v = variables.get(*i);

    // if not found, try to compute:
    if (v == NULL) {
      ierr = compute_by_name(*i, v); CHKERRQ(ierr);
    }

    // if still not found, ignore
    if (v == NULL)
      continue;

    GridType dims = v->grid_type();
    bool big_viewer = (big_viewers[*i]) || (big_viewers[*i + "_slice"]);

    // warn about 2D variables and ignore them:
    if (dims == GRID_2D) {
      ierr = verbPrintf(2, grid.com, "PISM WARNING: Please use -view instead of -view_slice to view 2D fields.\n");
      continue;
    }

    if (dims == GRID_3D) {
	IceModelVec3 *v3 = dynamic_cast<IceModelVec3*>(v);
	if (v3 == NULL) SETERRQ(1,"grid_type() gives GRID_3D but dynamic_cast gives a NULL");
	ierr = v3->view_horizontal_slice(slice_level, g2, big_viewer); CHKERRQ(ierr);
    }
  }

  // sounding viewers:
  // NOT IMPLEMENTED
  return 0;
}

PetscErrorCode IceModel::init_viewers() {
  PetscErrorCode ierr;
  PetscTruth flag;
  char tmp[TEMPORARY_STRING_LENGTH];

  ierr = PetscOptionsBegin(grid.com, PETSC_NULL,
			   "Options controlling run-time diagnostic viewers",
			   PETSC_NULL); CHKERRQ(ierr);

  // map-plane (and surface) viewers:
  ierr = PetscOptionsString("-view_map", "specifies the comma-separated list of map-plane viewers", "", "[empty]",
			    tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  string var_name;
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      map_viewers.insert(var_name);
  }

  // horizontal slice viewers:
  ierr = PetscOptionsString("-view_slice", "specifies the comma-separated list of horizontal-slice viewers", "", "empty",
			    tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      slice_viewers.insert(var_name);
  }

  // sounding viewers:
  ierr = PetscOptionsString("-view_sounding", "specifies the comma-separated list of sounding viewers", "", "empty",
			    tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      sounding_viewers.insert(var_name);
  }

  ierr = PetscOptionsString("-view_big", "specifies the comma-separated list of variables that should have big viewers",
			    "", "empty", tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      big_viewers[var_name] = true;
  }

  ierr = PetscOptionsReal("-slice_level", "sets the level (in meters above the base of ice) for slice viewers", "",
			  slice_level, &slice_level, PETSC_NULL); CHKERRQ(ierr);
  if ( (slice_level > grid.Lz) || (slice_level < 0) ) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: Slice level has to be positive and less than Lz (%3.3f).\n"
		      "              Disabling slice viewers...\n",
		      grid.Lz);
    slice_viewers.clear();
  }

  // Done with the options.
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

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
