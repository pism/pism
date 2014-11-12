// Copyright (C) 2004-2011, 2013, 2014 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include <petscdmda.h>
#include <petscksp.h>

#include "iceModel.hh"
#include "PISMDiagnostic.hh"

#include "error_handling.hh"

namespace pism {

//! Update the runtime graphical viewers.
/*!
Most viewers are updated by this routine, but some other are updated elsewhere.
 */
PetscErrorCode IceModel::update_viewers() {
  PetscErrorCode ierr;
  std::set<std::string>::iterator i;

  unsigned int viewer_size = (unsigned int)config.get("viewer_size");

  // map-plane viewers
  for (i = map_viewers.begin(); i != map_viewers.end(); ++i) {
    IceModelVec *v = variables.get(*i);
    bool de_allocate = false;

    // if not found, try to compute:
    if (v == NULL) {
      de_allocate = true;
      Diagnostic *diag = diagnostics[*i];

      if (diag) {
        ierr = diag->compute(v); CHKERRQ(ierr);
      } else {
        v = NULL;
      }
    }

    // if still not found, ignore
    if (v == NULL) {
      continue;
    }

    unsigned int dims = v->get_ndims();

    if (dims != 2) {
      throw RuntimeError("map-plane views of 3D quantities are not supported.");
    }

    if (v->get_ndof() == 1) {    // scalar fields
      std::string name = v->metadata().get_string("short_name");
      PetscViewer viewer = viewers[name];

      if (viewer == NULL) {
        ierr = grid.create_viewer(viewer_size, name, viewer); CHKERRQ(ierr);
        viewers[name] = viewer;
      }

      IceModelVec2S *v2d = dynamic_cast<IceModelVec2S*>(v);
      if (v2d == NULL) {
        throw RuntimeError("get_ndims() returns GRID_2D but dynamic_cast gives a NULL");
      }

      ierr = v2d->view(viewer, NULL); CHKERRQ(ierr);

    } else if (v->get_ndof() == 2) { // vector fields
      std::string name_1 = v->metadata().get_string("short_name"),
        name_2 = v->metadata(1).get_string("short_name");
      PetscViewer v1 = viewers[name_1],
        v2 = viewers[name_2];

      if (v1 == NULL) {
        ierr = grid.create_viewer(viewer_size, name_1, v1); CHKERRQ(ierr);
        viewers[name_1] = v1;
      }

      if (v2 == NULL) {
        ierr = grid.create_viewer(viewer_size, name_2, v2); CHKERRQ(ierr);
        viewers[name_2] = v2;
      }

      IceModelVec2 *v2d = dynamic_cast<IceModelVec2*>(v);
      if (v2d == NULL) {
        throw RuntimeError("get_ndims() returns GRID_2D but dynamic_cast gives a NULL");
      }

      ierr = v2d->view(v1, v2); CHKERRQ(ierr);
    }

    if (de_allocate) {
      delete v;
    }
  }

  return 0;
}

//! Initialize run-time diagnostic viewers.
PetscErrorCode IceModel::init_viewers() {
  PetscErrorCode ierr;
  PetscBool flag;
  char tmp[TEMPORARY_STRING_LENGTH];

  ierr = PetscOptionsBegin(grid.com, NULL,
                           "Options controlling run-time diagnostic viewers",
                           NULL);
  PISM_PETSC_CHK(ierr, "PetscOptionsBegin");

  PetscInt viewer_size = (int)config.get("viewer_size");
  ierr = PetscOptionsInt("-view_size", "specifies desired viewer size",
                         "", viewer_size, &viewer_size, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsInt");

  if (flag) {
    config.set_double("viewer_size", viewer_size);
  }

  // map-plane (and surface) viewers:
  ierr = PetscOptionsString("-view_map", "specifies the comma-separated list of map-plane viewers", "", "empty",
                            tmp, TEMPORARY_STRING_LENGTH, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");
  std::string var_name;
  if (flag) {
    std::istringstream arg(tmp);

    while (getline(arg, var_name, ',')) {
      map_viewers.insert(var_name);
    }
  }

  // Done with the options.
  ierr = PetscOptionsEnd();
  PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

  return 0;
}



} // end of namespace pism
