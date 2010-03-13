// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
Most viewers are updated by this routine, but some other are updated elsewhere:
  \li see computeMaxDiffusivity() in iMutil.cc for the diffusivity viewer.
  \li see update_nu_viewers() for nuH and log_nuH viewers.
 */
PetscErrorCode IceModel::update_viewers() {
  PetscErrorCode ierr;
  set<string>::iterator i;

  PetscInt viewer_size = (PetscInt)config.get("viewer_size");
  PetscScalar slice_level = config.get("slice_level");

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

    switch(dims) {
    case GRID_2D:
      {
	IceModelVec2 *v2 = dynamic_cast<IceModelVec2*>(v);
	if (v2 == NULL) SETERRQ(1,"grid_type() returns GRID_2D but dynamic_cast gives a NULL");
	ierr = v2->view(viewer_size); CHKERRQ(ierr);
	break;
      }
    case GRID_3D:
      {
	IceModelVec3 *v3 = dynamic_cast<IceModelVec3*>(v);
	if (v3 == NULL) SETERRQ(1,"grid_type() returns GRID_3D but dynamic_cast gives a NULL");
	ierr = v3->view_surface(vH, viewer_size); CHKERRQ(ierr);
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

    // warn about 2D variables and ignore them:
    if (dims == GRID_2D) {
      ierr = verbPrintf(2, grid.com, "PISM WARNING: Please use -view instead of -view_slice to view 2D fields.\n");
      continue;
    }

    if (dims == GRID_3D) {
	IceModelVec3 *v3 = dynamic_cast<IceModelVec3*>(v);
	if (v3 == NULL) SETERRQ(1,"grid_type() returns GRID_3D but dynamic_cast gives a NULL");
	ierr = v3->view_horizontal_slice(slice_level, viewer_size); CHKERRQ(ierr);
    }
  }

  // sounding viewers:
  for (i = sounding_viewers.begin(); i != sounding_viewers.end(); ++i) {
    IceModelVec *v = variables.get(*i);

    // if not found, try to compute:
    if (v == NULL) {
      ierr = compute_by_name(*i, v); CHKERRQ(ierr);
    }

    // if still not found, ignore
    if (v == NULL)
      continue;

    GridType dims = v->grid_type();

    // if it's a 2D variable, stop
    if (dims == GRID_2D) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: soundings of 2D quantities are not supported.\n");
      PetscEnd();
    }

    if (dims == GRID_3D) {
	IceModelVec3 *v3 = dynamic_cast<IceModelVec3*>(v);
	if (v3 == NULL) SETERRQ(1,"grid_type() returns GRID_3D but dynamic_cast gives a NULL");
	ierr = v3->view_sounding(id, jd, viewer_size); CHKERRQ(ierr);
    }

    if (dims == GRID_3D_BEDROCK) {
	IceModelVec3Bedrock *v3 = dynamic_cast<IceModelVec3Bedrock*>(v);
	if (v3 == NULL) SETERRQ(1,"grid_type() returns GRID_3D_BEDROCK but dynamic_cast gives a NULL");
	ierr = v3->view_sounding(id, jd, viewer_size); CHKERRQ(ierr);
    }
  }
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

  // map-plane (and surface) viewers:
  ierr = PetscOptionsString("-view_map", "specifies the comma-separated list of map-plane viewers", "", "empty",
			    tmp, TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);
  string var_name;
  if (flag) {
    istringstream arg(tmp);

    while (getline(arg, var_name, ',')) {
      if (var_name == "diffusivity")
	view_diffusivity = true;
      else if (var_name == "log_nuH")
	view_log_nuH = true;
      else if (var_name == "nuH")
	view_nuH = true;
      else
	map_viewers.insert(var_name);
    }
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

  PetscInt viewer_size = (PetscInt)config.get("viewer_size");
  ierr = PetscOptionsInt("-view_size", "specifies desired viewer size",
			 "", viewer_size, &viewer_size, &flag); CHKERRQ(ierr);

  if (flag)
    config.set("viewer_size", viewer_size); 

  PetscScalar slice_level = config.get("slice_level");
  ierr = PetscOptionsReal("-view_slice_level", "sets the level (in meters above the base of ice) for slice viewers", "",
			  slice_level, &slice_level, PETSC_NULL); CHKERRQ(ierr);
  if ( (slice_level > grid.Lz) || (slice_level < 0) ) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: Slice level has to be positive and less than Lz (%3.3f).\n"
		      "              Disabling slice viewers...\n",
		      grid.Lz);
    slice_viewers.clear();
  } else {
    config.set("slice_level", slice_level);
  }

  // Done with the options.
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  return 0;
}


//! Update nuH viewers.
PetscErrorCode IceModel::update_nu_viewers(IceModelVec2 vNu[2], IceModelVec2 /*vNuOld*/[2], bool /*updateNu_tView*/) {
  // this one is called when solving an SSA system

  PetscErrorCode ierr;
  if (view_log_nuH) {
    PetscScalar  **nui, **nuj, **gg;  
    ierr = vWork2d[5].get_array(gg); CHKERRQ(ierr);
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
    ierr = vWork2d[5].end_access(); CHKERRQ(ierr);

    ierr = vWork2d[5].set_name("log10(nuH)"); CHKERRQ(ierr);
    ierr = vWork2d[5].view(300); CHKERRQ(ierr);
  }

  if (view_nuH) {
    ierr = vNu[0].set_name("nuH[0]"); CHKERRQ(ierr);
    ierr = vNu[0].view(300); CHKERRQ(ierr);
    ierr = vNu[1].set_name("nuH[1]"); CHKERRQ(ierr);
    ierr = vNu[1].view(300); CHKERRQ(ierr);
  }

  return 0;
}

