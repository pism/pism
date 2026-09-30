// Copyright (C) 2026 PISM Authors
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

#include <utility> // std::move

#include "pism/coupler/surface/DEBMEnhanced.hh"
#include "pism/coupler/surface/DEBMSimplePointwise.hh"
#include "pism/coupler/surface/Insolation.hh"

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Context.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Array.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace surface {

///// dEBM-enhanced: dEBM-simple driven by a terrain-shaded surface-insolation field
///// computed from the ice surface elevation.

DEBMEnhanced::DEBMEnhanced(std::shared_ptr<const Grid> g,
                           std::shared_ptr<atmosphere::AtmosphereModel> input)
    : DEBMSimple(g, std::move(input)),
      m_update_interval(0.0),
      m_t_last_horizon(0.0) {

  m_log->message(2,
                 "  dEBM-enhanced: the terrain-shaded surface insolation is computed\n"
                 "  internally from the ice surface elevation.\n");

  DEBMSimpleAtmosphereTransmissivity transmissivity(*m_config);
  m_terrain.reset(new DEBMTerrainShadedInsolation(m_grid, transmissivity));

  m_update_interval =
      m_config->get_number("surface.debm_enhanced.update_interval", "seconds");
}

DEBMEnhanced::~DEBMEnhanced() = default;

void DEBMEnhanced::init_impl(const Inputs &inputs) {
  // initialize the dEBM-simple machinery first
  DEBMSimple::init_impl(inputs);

  m_log->message(2,
                 "* dEBM-enhanced: the insolation-driven melt uses terrain-shaded surface\n"
                 "  insolation computed from the ice surface elevation.\n");

  // Build the initial horizon map and sky view factor. The horizon is recomputed using
  // the modeled surface elevation every surface.debm_enhanced.update_interval.
  m_terrain->update_shading(inputs.geometry->ice_surface_elevation);
  m_t_last_horizon = m_grid->ctx()->time()->current();
}

void DEBMEnhanced::update_insolation_input(double t, double dt,
                                           const array::Scalar1 &surface_elevation,
                                           array::AccessScope &list) {

  // Recompute the terrain horizon from the current surface elevation once at least
  // m_update_interval has elapsed (a value of zero recomputes it every time step). The
  // horizon changes slowly as the geometry evolves, so this is much cheaper than redoing
  // the ray-marching every step.
  if (t >= m_t_last_horizon + m_update_interval) {
    m_terrain->update_shading(surface_elevation);
    m_t_last_horizon = t;
  }

  // Recompute the daily terrain-shaded insolation field for a representative day at the
  // midpoint of the update interval. Using the midpoint declination for the whole
  // interval is an approximation that is accurate for the sub-monthly time steps used in
  // practice.
  m_terrain->update_daily_insolation(t + 0.5 * dt, surface_elevation);
  list.add(m_terrain->insolation());
}

void DEBMEnhanced::insolation_energy_series(int i, int j,
                                            const std::vector<OrbitalParameters> &orbital,
                                            double dt_sub,
                                            double /*latitude*/,
                                            std::vector<double> &result) const {
  // the daily-mean insolation rate (W m^-2) at this cell was computed in
  // update_insolation_input
  double rate = m_terrain->insolation()(i, j);
  for (size_t k = 0; k < orbital.size(); ++k) {
    result[k] = rate * dt_sub;
  }
}

DiagnosticList DEBMEnhanced::spatial_diagnostics_impl() const {
  auto result = DEBMSimple::spatial_diagnostics_impl();
  // expose the insolation actually driving the melt, replacing the analytic dEBM-simple
  // "insolation" diagnostic (which dEBM-enhanced does not use)
  result["insolation"] = Diagnostic::wrap(m_terrain->insolation());
  // and the terrain horizon map (azimuth, y, x) and sky-view factor (y, x)
  result["horizon"] = Diagnostic::wrap(m_terrain->horizon_altitude());
  if (m_terrain->sky_view_enabled()) {
    result["sky_view_factor"] = Diagnostic::wrap(m_terrain->sky_view_factor());
  }
  return result;
}

} // end of namespace surface
} // end of namespace pism
