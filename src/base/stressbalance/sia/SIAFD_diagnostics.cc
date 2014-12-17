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
  dict["diffusivity"] = new SIAFD_diffusivity(this, m_grid);
  dict["diffusivity_staggered"] = new SIAFD_diffusivity_staggered(this, m_grid);
  dict["schoofs_theta"] = new SIAFD_schoofs_theta(this, m_grid);
  dict["thksmooth"] = new SIAFD_thksmooth(this, m_grid);
  dict["topgsmooth"] = new SIAFD_topgsmooth(this, m_grid);
  dict["h_x"] = new SIAFD_h_x(this, m_grid);
  dict["h_y"] = new SIAFD_h_y(this, m_grid);
}

SIAFD_schoofs_theta::SIAFD_schoofs_theta(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "schoofs_theta", grid));

  set_attrs("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA", "",
            "1", "", 0);
  vars[0].set_double("valid_min", 0);
  vars[0].set_double("valid_max", 1);
}

void SIAFD_schoofs_theta::compute(IceModelVec* &output) {
  IceModelVec2S *result, *surface;

  surface = grid.variables().get_2d_scalar("surface_altitude");

  result = new IceModelVec2S;
  result->create(grid, "schoofs_theta", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  model->bed_smoother->get_theta(*surface, result);

  output = result;
}


SIAFD_topgsmooth::SIAFD_topgsmooth(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "topgsmooth", grid));
  set_attrs("smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

void SIAFD_topgsmooth::compute(IceModelVec* &output) {
  IceModelVec2S *result;

  result = new IceModelVec2S;
  result->create(grid, "topgsmooth", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  result->copy_from(model->bed_smoother->get_smoothed_bed());

  output = result;
}

SIAFD_thksmooth::SIAFD_thksmooth(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "thksmooth", grid));
  set_attrs("thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

void SIAFD_thksmooth::compute(IceModelVec* &output) {
  IceModelVec2S *result, *surface, *thickness;
  IceModelVec2Int *mask;

  surface   = grid.variables().get_2d_scalar("surface_altitude");
  thickness = grid.variables().get_2d_scalar("land_ice_thickness");
  mask      = grid.variables().get_2d_mask("mask");

  result = new IceModelVec2S;
  result->create(grid, "thksmooth", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  model->bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                        result);

  output = result;
}



SIAFD_diffusivity::SIAFD_diffusivity(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "diffusivity", grid));

  set_attrs("diffusivity of SIA mass continuity equation", "",
            "m2 s-1", "m2 s-1", 0);
}

void SIAFD_diffusivity::compute(IceModelVec* &output) {
  IceModelVec2S *result;

  result = new IceModelVec2S;
  result->create(grid, "diffusivity", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  model->compute_diffusivity(*result);

  output = result;
}

SIAFD_diffusivity_staggered::SIAFD_diffusivity_staggered(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  dof = 2;

  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "diffusivity_i", grid));
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "diffusivity_j", grid));

  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (i-offset)", "",
            "m2 s-1", "m2 s-1", 0);
  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (j-offset)", "",
            "m2 s-1", "m2 s-1", 1);
}

void SIAFD_diffusivity_staggered::compute(IceModelVec* &output) {
  IceModelVec2Stag *result;

  result = new IceModelVec2Stag;
  result->create(grid, "diffusivity", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  model->compute_diffusivity_staggered(*result);

  output = result;
}

SIAFD_h_x::SIAFD_h_x(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  dof = 2;

  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "h_x_i", grid));
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "h_x_j", grid));

  set_attrs("the x-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the x-component of the surface gradient, j-offset", "",
            "", "", 1);
}

void SIAFD_h_x::compute(IceModelVec* &output) {

  IceModelVec2Stag *result = new IceModelVec2Stag;
  result->create(grid, "h_x", WITH_GHOSTS);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  model->compute_surface_gradient(model->work_2d_stag[0],
                                  model->work_2d_stag[1]);

  result->copy_from(model->work_2d_stag[0]);

  output = result;
}

SIAFD_h_y::SIAFD_h_y(SIAFD *m, IceGrid &g)
  : Diag<SIAFD>(m, g) {

  // set metadata:
  dof = 2;

  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "h_y_i", grid));
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "h_y_j", grid));

  set_attrs("the y-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the y-component of the surface gradient, j-offset", "",
            "", "", 1);
}

void SIAFD_h_y::compute(IceModelVec* &output) {

  IceModelVec2Stag *result = new IceModelVec2Stag;
  result->create(grid, "h_y", WITH_GHOSTS);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  model->compute_surface_gradient(model->work_2d_stag[0],
                                  model->work_2d_stag[1]);

  result->copy_from(model->work_2d_stag[1]);

  output = result;
}

} // end of namespace pism
