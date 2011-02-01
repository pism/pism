// Copyright (C) 2010, 2011 Constantine Khroulev
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

#include "SIAFD.hh"

void SIAFD::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  dict["diffusivity"] = new SIAFD_diffusivity(this, grid, *variables);
  dict["schoofs_theta"] = new SIAFD_schoofs_theta(this, grid, *variables);
  dict["thksmooth"] = new SIAFD_thksmooth(this, grid, *variables);
  dict["topgsmooth"] = new SIAFD_topgsmooth(this, grid, *variables);
}

SIAFD_schoofs_theta::SIAFD_schoofs_theta(SIAFD *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SIAFD>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init("schoofs_theta", grid, GRID_2D);
  
  set_attrs("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA", "",
            "1", "", 0);
  vars[0].set("valid_min", 0);
  vars[0].set("valid_max", 1);
}

PetscErrorCode SIAFD_schoofs_theta::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result, *surface;
  PetscInt WIDE_STENCIL = 2;

  result = new IceModelVec2S;
  ierr = result->create(grid, "schoofs_theta", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (surface == NULL) SETERRQ(1, "surface_altitude is not available");

  ierr = model->bed_smoother->get_theta(*surface, model->config.get("Glen_exponent"),
                                        WIDE_STENCIL, result); CHKERRQ(ierr);

  output = result;
  return 0;
}


SIAFD_topgsmooth::SIAFD_topgsmooth(SIAFD *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SIAFD>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init("topgsmooth", grid, GRID_2D);
  set_attrs("smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

PetscErrorCode SIAFD_topgsmooth::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result;
  PetscInt WIDE_STENCIL = 2;

  result = new IceModelVec2S;
  ierr = result->create(grid, "topgsmooth", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = result->copy_from(model->bed_smoother->topgsmooth); CHKERRQ(ierr);

  output = result;
  return 0;
}

SIAFD_thksmooth::SIAFD_thksmooth(SIAFD *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SIAFD>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init("thksmooth", grid, GRID_2D);
  set_attrs("thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

PetscErrorCode SIAFD_thksmooth::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = 2;
  IceModelVec2S *result, *surface, *thickness;
  IceModelVec2Mask *mask;

  result = new IceModelVec2S;
  ierr = result->create(grid, "thksmooth", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (surface == NULL) SETERRQ(1, "surface_altitude is not available");

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  mask = dynamic_cast<IceModelVec2Mask*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");

  ierr = model->bed_smoother->get_smoothed_thk(*surface, *thickness, *mask, WIDE_STENCIL,
                                               result); CHKERRQ(ierr);

  output = result;
  return 0;
}



SIAFD_diffusivity::SIAFD_diffusivity(SIAFD *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SIAFD>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init("diffusivity", grid, GRID_2D);
  
  set_attrs("diffusivity of SIA mass continuity equation", "",
            "m2 s-1", "m2 s-1", 0);
}

PetscErrorCode SIAFD_diffusivity::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result;
  
  result = new IceModelVec2S;
  ierr = result->create(grid, "diffusivity", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->compute_diffusivity(*result); CHKERRQ(ierr); 

  output = result;
  return 0;
}
