// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev
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

#include "SIAFD_diagnostics.hh"
#include "BedSmoother.hh"
#include "pism/util/Vars.hh"

namespace pism {
namespace stressbalance {

DiagnosticList SIAFD::diagnostics_impl() const {
  DiagnosticList result = {
    {"diffusivity",           Diagnostic::Ptr(new SIAFD_diffusivity(this))},
    {"diffusivity_staggered", Diagnostic::Ptr(new SIAFD_diffusivity_staggered(this))},
    {"schoofs_theta",         Diagnostic::Ptr(new SIAFD_schoofs_theta(this))},
    {"thksmooth",             Diagnostic::Ptr(new SIAFD_thksmooth(this))},
    {"topgsmooth",            Diagnostic::Ptr(new SIAFD_topgsmooth(this))},
    {"h_x",                   Diagnostic::Ptr(new SIAFD_h_x(this))},
    {"h_y",                   Diagnostic::Ptr(new SIAFD_h_y(this))}
  };
  return result;
}

SIAFD_schoofs_theta::SIAFD_schoofs_theta(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "schoofs_theta")};

  set_attrs("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA", "",
            "1", "", 0);
  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("valid_max", 1);
}

IceModelVec::Ptr SIAFD_schoofs_theta::compute_impl() const {
  const IceModelVec2S *surface = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "schoofs_theta", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->bed_smoother().theta(*surface, *result);

  return result;
}


SIAFD_topgsmooth::SIAFD_topgsmooth(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "topgsmooth")};
  set_attrs("smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

IceModelVec::Ptr SIAFD_topgsmooth::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "topgsmooth", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  result->copy_from(model->bed_smoother().smoothed_bed());

  return result;
}

SIAFD_thksmooth::SIAFD_thksmooth(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "thksmooth")};
  set_attrs("thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

IceModelVec::Ptr SIAFD_thksmooth::compute_impl() const {

  const IceModelVec2S        &surface   = *m_grid->variables().get_2d_scalar("surface_altitude");
  const IceModelVec2S        &thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2CellType &mask      = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "thksmooth", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->bed_smoother().smoothed_thk(surface, thickness, mask,
                                     *result);
  return result;
}



SIAFD_diffusivity::SIAFD_diffusivity(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "diffusivity")};

  set_attrs("diffusivity of SIA mass continuity equation", "",
            "m2 s-1", "m2 s-1", 0);
}

IceModelVec::Ptr SIAFD_diffusivity::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "diffusivity", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  model->diffusivity().staggered_to_regular(*result);

  return result;
}

SIAFD_diffusivity_staggered::SIAFD_diffusivity_staggered(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "diffusivity_i"),
            SpatialVariableMetadata(m_sys, "diffusivity_j")};

  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (i-offset)", "",
            "m2 s-1", "m2 s-1", 0);
  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (j-offset)", "",
            "m2 s-1", "m2 s-1", 1);
}

static void copy_staggered_vec(const IceModelVec2Stag &input, IceModelVec2Stag &output) {
  IceGrid::ConstPtr grid = output.grid();

  IceModelVec::AccessList list{ &input, &output };

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    output(i, j, 0) = input(i, j, 0);
    output(i, j, 1) = input(i, j, 1);
  }
}

IceModelVec::Ptr SIAFD_diffusivity_staggered::compute_impl() const {
  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "diffusivity", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  copy_staggered_vec(model->diffusivity(), *result.get());

  return result;
}

SIAFD_h_x::SIAFD_h_x(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "h_x_i"),
            SpatialVariableMetadata(m_sys, "h_x_j")};

  set_attrs("the x-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the x-component of the surface gradient, j-offset", "",
            "", "", 1);
}

IceModelVec::Ptr SIAFD_h_x::compute_impl() const {

  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "h_x", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  copy_staggered_vec(model->surface_gradient_x(), *result.get());

  return result;
}

SIAFD_h_y::SIAFD_h_y(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "h_y_i"),
            SpatialVariableMetadata(m_sys, "h_y_j")};

  set_attrs("the y-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the y-component of the surface gradient, j-offset", "",
            "", "", 1);
}

IceModelVec::Ptr SIAFD_h_y::compute_impl() const {

  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "h_y", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  copy_staggered_vec(model->surface_gradient_y(), *result.get());

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
