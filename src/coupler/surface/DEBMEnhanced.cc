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
#include "pism/coupler/surface/TerrainInsolation.hh"

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
#include "pism/util/error_handling.hh"

namespace pism {
namespace surface {

///// dEBM-enhanced: dEBM-simple driven by a terrain-shaded surface-insolation field
///// computed from the ice surface elevation.

DEBMEnhanced::DEBMEnhanced(std::shared_ptr<const Grid> g,
                           std::shared_ptr<atmosphere::AtmosphereModel> input)
  : DEBMSimple(g, std::move(input)), m_latitude(nullptr), m_surface_elevation(nullptr),
    m_update_interval(0.0), m_t_last_horizon(0.0) {

  // dEBM-enhanced reuses the analytic orbit only for the temperature/offset melt-period
  // weighting. The insolation is computed for a fixed (present-day) orbit, so it is
  // inconsistent with the paleo orbit parameterization.
  if (m_config->get_flag("surface.debm_simple.paleo.enabled")) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "the 'debm_enhanced' surface model is incompatible with "
                                  "surface.debm_simple.paleo.enabled");
  }

  m_log->message(2,
                 "  dEBM-enhanced: the terrain-shaded surface insolation is computed\n"
                 "  internally from the ice surface elevation.\n");

  m_terrain.reset(new TerrainInsolation(m_grid));

  m_update_interval =
      m_config->get_number("surface.debm_enhanced.update_interval", "seconds");

  m_computed_insolation = std::make_shared<array::Scalar>(m_grid, "insolation");
  m_computed_insolation->metadata(0)
      .long_name("daily mean terrain-shaded surface insolation")
      .units("W m^-2");
}

// Defaulted here (not in the header) because m_terrain is a unique_ptr to the
// forward-declared TerrainInsolation, whose complete type is needed by the deleter.
DEBMEnhanced::~DEBMEnhanced() = default;

void DEBMEnhanced::init_impl(const Geometry &geometry) {
  // initialize the dEBM-simple machinery first
  DEBMSimple::init_impl(geometry);

  m_log->message(2,
                 "* dEBM-enhanced: the insolation-driven melt uses terrain-shaded surface\n"
                 "  insolation computed from the ice surface elevation (instead of the\n"
                 "  analytic top-of-atmosphere parameterization).\n");

  // Cache the latitude and surface-elevation fields and build the initial horizon map and
  // surface normals. The horizon is recomputed from the (evolving) surface elevation every
  // surface.debm_enhanced.update_interval (see update_insolation_input).
  m_latitude = &geometry.latitude;
  m_surface_elevation = &geometry.ice_surface_elevation;
  m_terrain->init(geometry.ice_surface_elevation);
  m_t_last_horizon = m_grid->ctx()->time()->current();
}

void DEBMEnhanced::update_insolation_input(double t, double dt,
                                           const std::vector<double> &ts,
                                           array::AccessScope &list) {
  (void)ts;

  // Recompute the terrain horizon from the current ice surface elevation once at least
  // m_update_interval has elapsed (a value of zero recomputes it every time step). The
  // horizon changes slowly as the geometry evolves, so this is much cheaper than redoing
  // the ray-marching every step.
  if (t >= m_t_last_horizon + m_update_interval) {
    m_terrain->init(*m_surface_elevation);
    m_t_last_horizon = t;
  }

  // Recompute the daily terrain-shaded insolation field for a representative day at the
  // midpoint of the update interval. The diurnal integral is independent of longitude, so
  // only the day's solar declination and Sun-Earth distance factor are needed (plus the
  // static horizon and normals). Using the midpoint declination for the whole interval is
  // an approximation that is accurate for the sub-monthly time steps used in practice.
  double t_mid = t + 0.5 * dt;
  double year_fraction = m_grid->ctx()->time()->year_fraction(t_mid);
  double declination = DEBMSimplePointwise::solar_declination_present_day(year_fraction);
  double distance_factor = DEBMSimplePointwise::distance_factor_present_day(year_fraction);

  m_terrain->daily_insolation(declination, distance_factor, *m_latitude,
                              *m_computed_insolation);
  list.add(*m_computed_insolation);
}

void DEBMEnhanced::insolation_energy_series(int i, int j,
                                            const std::vector<DEBMSimpleOrbitalParameters> &orbital,
                                            const std::vector<double> &ts,
                                            double dt_sub,
                                            double latitude,
                                            std::vector<double> &result) const {
  (void)ts;       // sampling times are already set via init_interpolation()
  (void)latitude; // spatial variability is encoded in the computed field

  // the daily-mean insolation rate (W m^-2) at this cell was computed in
  // update_insolation_input; the energy reaching the surface during a sub-step of length
  // dt_sub (seconds) is rate * dt_sub
  double rate = (*m_computed_insolation)(i, j);
  for (size_t k = 0; k < orbital.size(); ++k) {
    result[k] = rate * dt_sub;
  }
}

DiagnosticList DEBMEnhanced::spatial_diagnostics_impl() const {
  DiagnosticList result = DEBMSimple::spatial_diagnostics_impl();
  // expose the insolation actually driving the melt, replacing the analytic dEBM-simple
  // "insolation" diagnostic (which dEBM-enhanced does not use)
  result["insolation"] = Diagnostic::wrap(*m_computed_insolation);
  // and the terrain horizon map (azimuth, y, x) and sky-view factor (y, x)
  result["horizon"] = Diagnostic::wrap(m_terrain->horizon());
  if (m_terrain->sky_view_enabled()) {
    result["sky_view_factor"] = Diagnostic::wrap(m_terrain->sky_view());
  }
  return result;
}

} // end of namespace surface
} // end of namespace pism
