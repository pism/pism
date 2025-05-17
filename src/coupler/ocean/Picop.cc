// Copyright (C) 2012-2019, 2021, 2022, 2023, 2024, 2025 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, and Matthias Mengel
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
//
// Please cite this model as:
// 1.
// Antarctic sub-shelf melt rates via PICO
// R. Reese, T. Albrecht, M. Mengel, X. Asay-Davis and R. Winkelmann
// The Cryosphere, 12, 1969-1985, (2018)
// DOI: 10.5194/tc-12-1969-2018
//
// 2.
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z
//
// 3.
// PICOP, a new ocean melt parameterization under ice shelves
// combining PICO and a plume model.
// T. Pelle, M. Morlighem, J.H. Bondzio
// The Cryosphere, 13, 1043-49, (2019)
// DOI: 10.5194/tc-13-1043-2019

#include <gsl/gsl_math.h> // GSL_NAN

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"

#include "pism/coupler/ocean/Picop.hh"
#include "pism/coupler/ocean/PicoGeometry.hh"
#include "pism/coupler/ocean/PicoPhysics.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace ocean {

Picop::Picop(std::shared_ptr<const Grid> grid)
  : CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),
    m_pico(std::make_shared<Pico>(grid)),
    m_grounding_line_elevation(grid, "picop_grounding_line_elevation"),
    m_geometry(grid),
    m_velocity(grid, "ghosted_velocity")
{

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  {
    auto buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_theta_ocean = std::make_shared<array::Forcing>(m_grid,
                                                file,
                                                "theta_ocean",
                                                "", // no standard name
                                                buffer_size,
                                                opt.periodic,
                                                LINEAR);

    m_salinity_ocean = std::make_shared<array::Forcing>(m_grid,
                                                   file,
                                                   "salinity_ocean",
                                                   "", // no standard name
                                                   buffer_size,
                                                   opt.periodic,
                                                   LINEAR);
  }

  m_theta_ocean->metadata(0)
      .long_name("potential temperature of the adjacent ocean")
      .units("kelvin");

  m_salinity_ocean->metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg");


  m_grounding_line_elevation.metadata(0).long_name("cavity overturning").units("m");
  m_grounding_line_elevation.metadata()["_FillValue"] = { 0.0 };
  
}

void Picop::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the Potsdam Ice-shelf Cavity mOdel / Plume for the ocean ...\n");

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  m_theta_ocean->init(opt.filename, opt.periodic);
  m_salinity_ocean->init(opt.filename, opt.periodic);

  // This initializes the basin_mask
  m_geometry.init();

  compute_grounding_line_elevation(geometry, m_grounding_line_elevation);
}

void Picop::define_model_state_impl(const File &output) const {


  OceanModel::define_model_state_impl(output);
}

void Picop::write_model_state_impl(const File &output) const {


  OceanModel::write_model_state_impl(output);
}


void Picop::update_impl(const Geometry &geometry, double t, double dt) {

  m_theta_ocean->update(t, dt);
  m_salinity_ocean->update(t, dt);

  m_theta_ocean->average(t, dt);
  m_salinity_ocean->average(t, dt);

  compute_grounding_line_elevation(geometry, m_grounding_line_elevation);
}


MaxTimestep Picop::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("ocean picop");
}


void Picop::compute_grounding_line_elevation(const Geometry &geometry,
                                             array::Scalar &grounding_line_elevation) const {
  const array::Scalar &bed      = geometry.bed_elevation;
  const array::Scalar &H        = geometry.ice_thickness;
  const array::Scalar &cell_type = geometry.cell_type;
  const array::Scalar &z_s      = geometry.sea_level_elevation;

  array::AccessScope list{&bed, &H, &z_s, &cell_type, &m_velocity, &grounding_line_elevation};

  auto grid = m_grid;

  // Step 1: Initialize zgl0 at grounding line: bed elevation
  for (auto p = grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    if (cell_type(i, j) == MASK_GROUNDED) {
      grounding_line_elevation(i, j) = bed(i, j);
    } else {
      grounding_line_elevation(i, j) = 0.0; // temporary initial value
    }
  }

  // Step 2: Advection-diffusion loop (simplified iterative update)
  const double eps = 1e-14;
  const int max_iter = 100;
  const double tol = 1e-4;
  double residual = 0.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    residual = 0.0;
    for (auto p = grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();
      if (cell_type(i, j) != MASK_FLOATING) continue;

      double u = m_velocity(i, j).u;
      double v = m_velocity(i, j).v;

      // Upwind scheme for advection term
      double dzdx = (u >= 0.0) ? (grounding_line_elevation(i, j) - grounding_line_elevation(i - 1, j)) / grid->dx()
                                : (grounding_line_elevation(i + 1, j) - grounding_line_elevation(i, j)) / grid->dx();
      double dzdy = (v >= 0.0) ? (grounding_line_elevation(i, j) - grounding_line_elevation(i, j - 1)) / grid->dy()
                                : (grounding_line_elevation(i, j + 1) - grounding_line_elevation(i, j)) / grid->dy();

      // Diffusion (Laplacian)
      double laplacian =
          (grounding_line_elevation(i + 1, j) + grounding_line_elevation(i - 1, j)
         + grounding_line_elevation(i, j + 1) + grounding_line_elevation(i, j - 1)
         - 4.0 * grounding_line_elevation(i, j)) / (grid->dx() * grid->dx());

      // Update
      double update = - (u * dzdx + v * dzdy + eps * laplacian);
      grounding_line_elevation(i, j) += 0.1 * update;

      residual += update * update;
    }

    residual = std::sqrt(GlobalSum(grid->com, residual));
    if (residual < tol) break;
  }

  // Step 3: Clip to ensure z_gl <= base of ice shelf
  for (auto p = grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    double zb = z_s(i, j) - H(i, j); // ice shelf base
    if (grounding_line_elevation(i, j) > zb) {
      grounding_line_elevation(i, j) = zb;
    }
  }

  m_log->message(2, "Computed grounding line elevation field with residual = %.3e\n", residual);
}

// Write diagnostic variables to extra files if requested
DiagnosticList Picop::diagnostics_impl() const {

  DiagnosticList result = {
    { "picop_grounding_line_elevation", Diagnostic::wrap(m_grounding_line_elevation) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}


} // end of namespace ocean
} // end of namespace pism
