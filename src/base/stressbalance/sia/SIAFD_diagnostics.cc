// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
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

#include "SIAFD.hh"
#include "PISMBedSmoother.hh"
#include "PISMVars.hh"

namespace pism {

void SIAFD::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                            std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  dict["diffusivity"] = new SIAFD_diffusivity(this, grid, *variables);
  dict["diffusivity_staggered"] = new SIAFD_diffusivity_staggered(this, grid, *variables);
  dict["schoofs_theta"] = new SIAFD_schoofs_theta(this, grid, *variables);
  dict["thksmooth"] = new SIAFD_thksmooth(this, grid, *variables);
  dict["topgsmooth"] = new SIAFD_topgsmooth(this, grid, *variables);
  dict["h_x"] = new SIAFD_h_x(this, grid, *variables);
  dict["h_y"] = new SIAFD_h_y(this, grid, *variables);
}

SIAFD_schoofs_theta::SIAFD_schoofs_theta(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("schoofs_theta", grid);

  set_attrs("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA", "",
            "1", "", 0);
  vars[0].set_double("valid_min", 0);
  vars[0].set_double("valid_max", 1);
}

PetscErrorCode SIAFD_schoofs_theta::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result, *surface;

  surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (surface == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  result = new IceModelVec2S;
  ierr = result->create(grid, "schoofs_theta", WITH_GHOSTS,
                        surface->get_stencil_width()); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->bed_smoother->get_theta(*surface, result); CHKERRQ(ierr);

  output = result;
  return 0;
}


SIAFD_topgsmooth::SIAFD_topgsmooth(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("topgsmooth", grid);
  set_attrs("smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

PetscErrorCode SIAFD_topgsmooth::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result;

  result = new IceModelVec2S;
  ierr = result->create(grid, "topgsmooth", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = result->copy_from(model->bed_smoother->get_smoothed_bed()); CHKERRQ(ierr);

  output = result;
  return 0;
}

SIAFD_thksmooth::SIAFD_thksmooth(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("thksmooth", grid);
  set_attrs("thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

PetscErrorCode SIAFD_thksmooth::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result, *surface, *thickness;
  IceModelVec2Int *mask;

  surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (surface == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  result = new IceModelVec2S;
  ierr = result->create(grid, "thksmooth", WITH_GHOSTS,
                        surface->get_stencil_width()); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                               result); CHKERRQ(ierr);

  output = result;
  return 0;
}



SIAFD_diffusivity::SIAFD_diffusivity(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("diffusivity", grid);

  set_attrs("diffusivity of SIA mass continuity equation", "",
            "m2 s-1", "m2 s-1", 0);
}

PetscErrorCode SIAFD_diffusivity::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result;

  result = new IceModelVec2S;
  ierr = result->create(grid, "diffusivity", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->compute_diffusivity(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

SIAFD_diffusivity_staggered::SIAFD_diffusivity_staggered(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof, NCSpatialVariable(g.get_unit_system()));
  vars[0].init_2d("diffusivity_i", grid);
  vars[1].init_2d("diffusivity_j", grid);

  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (i-offset)", "",
            "m2 s-1", "m2 s-1", 0);
  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (j-offset)", "",
            "m2 s-1", "m2 s-1", 1);
}

PetscErrorCode SIAFD_diffusivity_staggered::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2Stag *result;

  result = new IceModelVec2Stag;
  ierr = result->create(grid, "diffusivity", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  ierr = model->compute_diffusivity_staggered(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

SIAFD_h_x::SIAFD_h_x(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof, NCSpatialVariable(g.get_unit_system()));
  vars[0].init_2d("h_x_i", grid);
  vars[1].init_2d("h_x_j", grid);

  set_attrs("the x-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the x-component of the surface gradient, j-offset", "",
            "", "", 1);
}

PetscErrorCode SIAFD_h_x::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2Stag *result = new IceModelVec2Stag;
  ierr = result->create(grid, "h_x", WITH_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  ierr = model->compute_surface_gradient(model->work_2d_stag[0],
                                         model->work_2d_stag[1]); CHKERRQ(ierr);

  ierr = result->copy_from(model->work_2d_stag[0]); CHKERRQ(ierr);

  output = result;
  return 0;
}

SIAFD_h_y::SIAFD_h_y(SIAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SIAFD>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof, NCSpatialVariable(g.get_unit_system()));
  vars[0].init_2d("h_y_i", grid);
  vars[1].init_2d("h_y_j", grid);

  set_attrs("the y-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the y-component of the surface gradient, j-offset", "",
            "", "", 1);
}

PetscErrorCode SIAFD_h_y::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2Stag *result = new IceModelVec2Stag;
  ierr = result->create(grid, "h_y", WITH_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  ierr = model->compute_surface_gradient(model->work_2d_stag[0],
                                         model->work_2d_stag[1]); CHKERRQ(ierr);

  ierr = result->copy_from(model->work_2d_stag[1]); CHKERRQ(ierr);

  output = result;
  return 0;
}

} // end of namespace pism
