// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

//! Update the runtime graphical viewers.
/*!
Most viewers are updated by this routine, but some other are updated elsewhere.
 */
PetscErrorCode IceModel::update_viewers() {
  PetscErrorCode ierr;
  set<string>::iterator i;

  PetscInt viewer_size = (PetscInt)config.get("viewer_size");

  // map-plane viewers
  for (i = map_viewers.begin(); i != map_viewers.end(); ++i) {
    IceModelVec *v = variables.get(*i);
    bool de_allocate = false;

    // if not found, try to compute:
    if (v == NULL) {
      de_allocate = true;
      PISMDiagnostic *diag = diagnostics[*i];

      if (diag) {
        ierr = diag->compute(v); CHKERRQ(ierr);
      } else {
        v = NULL;
      }
    }

    // if still not found, ignore
    if (v == NULL)
      continue;

    GridType dims = v->grid_type();

    if (dims != GRID_2D) {
      ierr = PetscPrintf(grid.com,
                         "PISM ERROR: map-plane views of 3D quantities are not supported.\n");
      CHKERRQ(ierr);
      PISMEnd();
    }

    if (v->get_dof() == 1) {    // scalar fields
      string name = v->string_attr("short_name");
      PetscViewer viewer = viewers[name];

      if (viewer == PETSC_NULL) {
        ierr = grid.create_viewer(viewer_size, name, viewer); CHKERRQ(ierr);
        viewers[name] = viewer;
      }

      IceModelVec2S *v2d = dynamic_cast<IceModelVec2S*>(v);
      if (v2d == NULL) SETERRQ(1,"grid_type() returns GRID_2D but dynamic_cast gives a NULL");

      ierr = v2d->view(viewer); CHKERRQ(ierr);

    } else if (v->get_dof() == 2) { // vector fields
      string name_1 = v->string_attr("short_name", 0),
        name_2 = v->string_attr("short_name", 1);
      PetscViewer v1 = viewers[name_1],
        v2 = viewers[name_2];

      if (v1 == PETSC_NULL) {
        ierr = grid.create_viewer(viewer_size, name_1, v1); CHKERRQ(ierr);
        viewers[name_1] = v1;
      }

      if (v2 == PETSC_NULL) {
        ierr = grid.create_viewer(viewer_size, name_2, v2); CHKERRQ(ierr);
        viewers[name_2] = v2;
      }

      IceModelVec2V *v2d = dynamic_cast<IceModelVec2V*>(v);
      if (v2d == NULL) SETERRQ(1,"grid_type() returns GRID_2D but dynamic_cast gives a NULL");

      ierr = v2d->view(v1, v2); CHKERRQ(ierr);
    }
    
    if (de_allocate) delete v;
  }

  // sounding viewers:
  for (i = sounding_viewers.begin(); i != sounding_viewers.end(); ++i) {
    IceModelVec *v = variables.get(*i);
    bool de_allocate = false;

    // if not found, try to compute:
    if (v == NULL) {
      de_allocate = true;
      PISMDiagnostic *diag = diagnostics[*i];

      if (diag) {
        ierr = diag->compute(v); CHKERRQ(ierr);
      } else {
        v = NULL;
      }
    }

    // if still not found, ignore
    if (v == NULL)
      continue;

    GridType dims = v->grid_type();

    // if it's a 2D variable, stop
    if (dims == GRID_2D) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: soundings of 2D quantities are not supported.\n");
      PISMEnd();
    }

    string name = v->string_attr("short_name");
    PetscViewer viewer = viewers[name];

    if (viewer == PETSC_NULL) {
      ierr = grid.create_viewer(viewer_size, name, viewers[name]); CHKERRQ(ierr);
      viewer = viewers[name];
    }

    if (dims == GRID_3D) {
	IceModelVec3 *v3d = dynamic_cast<IceModelVec3*>(v);
	if (v3d == NULL) SETERRQ(1,"grid_type() returns GRID_3D but dynamic_cast gives a NULL");
	ierr = v3d->view_sounding(id, jd, viewer); CHKERRQ(ierr);
    }

    if (dims == GRID_3D_BEDROCK) {
	IceModelVec3Bedrock *v3d = dynamic_cast<IceModelVec3Bedrock*>(v);
	if (v3d == NULL) SETERRQ(1,"grid_type() returns GRID_3D_BEDROCK but dynamic_cast gives a NULL");
	ierr = v3d->view_sounding(id, jd, viewer); CHKERRQ(ierr);
    }

    if (de_allocate) delete v;
  } // sounding viewers
  return 0;
}

//! Initialize run-time diagnostic viewers.
PetscErrorCode IceModel::init_viewers() {
  PetscErrorCode ierr;
  PetscTruth flag;
  char tmp[TEMPORARY_STRING_LENGTH];

  ierr = PetscOptionsBegin(grid.com, PETSC_NULL,
			   "Options controlling run-time diagnostic viewers",
			   PETSC_NULL); CHKERRQ(ierr);

  PetscInt viewer_size = (PetscInt)config.get("viewer_size");
  ierr = PetscOptionsInt("-view_size", "specifies desired viewer size",
			 "", viewer_size, &viewer_size, &flag); CHKERRQ(ierr);

  if (flag)
    config.set("viewer_size", viewer_size); 

  // map-plane (and surface) viewers:
  ierr = PetscOptionsString("-view_map", "specifies the comma-separated list of map-plane viewers", "", "empty",
			    tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  string var_name;
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ',')) {
      map_viewers.insert(var_name);
    }
  }

  // sounding viewers:
  ierr = PetscOptionsString("-view_sounding", "specifies the comma-separated list of sounding viewers", "", "empty",
			    tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ','))
      sounding_viewers.insert(var_name);
  }

  // Done with the options.
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  return 0;
}


