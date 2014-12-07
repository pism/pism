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
void IceModel::update_viewers() {
  std::set<std::string>::iterator i;

  unsigned int viewer_size = (unsigned int)config.get("viewer_size");

  // map-plane viewers
  for (i = map_viewers.begin(); i != map_viewers.end(); ++i) {
    IceModelVec *v = NULL;
    bool de_allocate = false;

    // if not found, try to compute:
    if (not variables.is_available(*i)) {
      de_allocate = true;
      Diagnostic *diag = diagnostics[*i];

      if (diag) {
        diag->compute(v);
      } else {
        v = NULL;
      }
    } else {
      v = variables.get(*i);
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
      Viewer::Ptr viewer = viewers[name];

      if (not viewer) {
        viewers[name].reset(new Viewer(grid.com, name, viewer_size, grid.Lx(), grid.Ly()));
        viewer = viewers[name];
      }

      IceModelVec2S *v2d = dynamic_cast<IceModelVec2S*>(v);
      if (v2d == NULL) {
        throw RuntimeError("get_ndims() returns GRID_2D but dynamic_cast gives a NULL");
      }

      v2d->view(viewer, Viewer::Ptr());

    } else if (v->get_ndof() == 2) { // vector fields
      std::string
        name_1 = v->metadata(0).get_string("short_name"),
        name_2 = v->metadata(1).get_string("short_name");
      Viewer::Ptr
        v1 = viewers[name_1],
        v2 = viewers[name_2];

      if (not v1) {
        viewers[name_1].reset(new Viewer(grid.com, name_1, viewer_size, grid.Lx(), grid.Ly()));
        v1 = viewers[name_1];
      }

      if (not v2) {
        viewers[name_2].reset(new Viewer(grid.com, name_2, viewer_size, grid.Lx(), grid.Ly()));
        v2 = viewers[name_2];
      }

      IceModelVec2 *v2d = dynamic_cast<IceModelVec2*>(v);
      if (v2d == NULL) {
        throw RuntimeError("get_ndims() returns GRID_2D but dynamic_cast gives a NULL");
      }

      v2d->view(v1, v2);
    }

    if (de_allocate) {
      delete v;
    }
  }
}

//! Initialize run-time diagnostic viewers.
void IceModel::init_viewers() {
  PetscErrorCode ierr;
  PetscBool flag;
  char tmp[TEMPORARY_STRING_LENGTH];

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
}



} // end of namespace pism
