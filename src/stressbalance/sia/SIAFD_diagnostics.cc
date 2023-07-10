// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2021, 2022, 2023 Constantine Khroulev
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

#include "pism/stressbalance/sia/SIAFD_diagnostics.hh"
#include "pism/stressbalance/sia/BedSmoother.hh"
#include "pism/util/Vars.hh"
#include "pism/util/array/CellType.hh"

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
  m_vars = { { m_sys, "schoofs_theta" } };

  m_vars[0]
      .long_name("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA")
      .units("1");
  m_vars[0]["valid_range"] = {0.0, 1.0};
}

std::shared_ptr<array::Array> SIAFD_schoofs_theta::compute_impl() const {
  const array::Scalar *surface = m_grid->variables().get_2d_scalar("surface_altitude");

  std::shared_ptr<array::Scalar> result(new array::Scalar(m_grid, "schoofs_theta"));
  result->metadata(0) = m_vars[0];

  model->bed_smoother().theta(*surface, *result);

  return result;
}


SIAFD_topgsmooth::SIAFD_topgsmooth(const SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars = { { m_sys, "topgsmooth" } };
  m_vars[0]
      .long_name("smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA")
      .units("m");
}

std::shared_ptr<array::Array> SIAFD_topgsmooth::compute_impl() const {

  std::shared_ptr<array::Scalar> result(new array::Scalar(m_grid, "topgsmooth"));
  result->metadata() = m_vars[0];

  result->copy_from(model->bed_smoother().smoothed_bed());

  return result;
}

SIAFD_thksmooth::SIAFD_thksmooth(const SIAFD *m)
  : Diag<SIAFD>(m) {

  m_vars = { { m_sys, "thksmooth" } };
  m_vars[0]
      .long_name(
          "thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA")
      .units("m");
}

std::shared_ptr<array::Array> SIAFD_thksmooth::compute_impl() const {

  const auto &surface   = *m_grid->variables().get_2d_scalar("surface_altitude");
  const auto &thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  array::CellType2 cell_type(m_grid, "cell_type");
  {
    const auto &mask = *m_grid->variables().get_2d_cell_type("mask");
    cell_type.copy_from(mask);
  }

  std::shared_ptr<array::Scalar> result(new array::Scalar(m_grid, "thksmooth"));
  result->metadata(0) = m_vars[0];

  model->bed_smoother().smoothed_thk(surface, thickness, cell_type,
                                     *result);
  return result;
}



SIAFD_diffusivity::SIAFD_diffusivity(const SIAFD *m)
  : Diag<SIAFD>(m) {

  m_vars = { { m_sys, "diffusivity" } };
  m_vars[0].long_name("diffusivity of SIA mass continuity equation").units("m2 s-1");
}

std::shared_ptr<array::Array> SIAFD_diffusivity::compute_impl() const {
  std::shared_ptr<array::Scalar> result(new array::Scalar(m_grid, "diffusivity"));
  result->metadata() = m_vars[0];

  array::CellType1 cell_type(m_grid, "cell_type");
  {
    const auto &mask = *m_grid->variables().get_2d_cell_type("mask");
    cell_type.copy_from(mask);
  }
  bool include_floating_ice = true;
  staggered_to_regular(cell_type, model->diffusivity(), include_floating_ice, *result);

  return result;
}

SIAFD_diffusivity_staggered::SIAFD_diffusivity_staggered(const SIAFD *m)
  : Diag<SIAFD>(m) {

  m_vars = { { m_sys, "diffusivity_i" }, { m_sys, "diffusivity_j" } };
  m_vars[0]
      .long_name("diffusivity of SIA mass continuity equation on the staggered grid (i-offset)")
      .units("m2 s-1");
  m_vars[1]
      .long_name("diffusivity of SIA mass continuity equation on the staggered grid (j-offset)")
      .units("m2 s-1");
}

static void copy_staggered_vec(const array::Staggered &input, array::Staggered &output) {
  auto grid = output.grid();

  array::AccessScope list{ &input, &output };

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    output(i, j, 0) = input(i, j, 0);
    output(i, j, 1) = input(i, j, 1);
  }
}

std::shared_ptr<array::Array> SIAFD_diffusivity_staggered::compute_impl() const {
  auto result = std::make_shared<array::Staggered>(m_grid, "diffusivity");

  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  copy_staggered_vec(model->diffusivity(), *result);

  return result;
}

SIAFD_h_x::SIAFD_h_x(const SIAFD *m)
  : Diag<SIAFD>(m) {

  m_vars = { { m_sys, "h_x_i" }, { m_sys, "h_x_j" } };
  m_vars[0].long_name("the x-component of the surface gradient, i-offset").units("1");
  m_vars[1].long_name("the x-component of the surface gradient, j-offset").units("1");
}

std::shared_ptr<array::Array> SIAFD_h_x::compute_impl() const {

  auto result = std::make_shared<array::Staggered>(m_grid, "h_x");

  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  copy_staggered_vec(model->surface_gradient_x(), *result.get());

  return result;
}

SIAFD_h_y::SIAFD_h_y(const SIAFD *m)
  : Diag<SIAFD>(m) {

  m_vars = { { m_sys, "h_y_i" }, { m_sys, "h_y_j" } };
  m_vars[0].long_name("the y-component of the surface gradient, i-offset").units("1");
  m_vars[1].long_name("the y-component of the surface gradient, j-offset").units("1");
}

std::shared_ptr<array::Array> SIAFD_h_y::compute_impl() const {

  auto result = std::make_shared<array::Staggered>(m_grid, "h_y");

  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  copy_staggered_vec(model->surface_gradient_y(), *result.get());

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
