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

#include "pism/coupler/atmosphere/CARRA2Bias.hh"

#include <algorithm> // std::min, std::max
#include <cmath>     // std::exp, std::cos, std::sin, std::ceil

#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace atmosphere {

CARRA2Bias::CARRA2Bias(std::shared_ptr<const Grid> grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in),
    m_surface(grid, "ice_surface_elevation") {

  const std::string prefix = "atmosphere.carra2_bias.";

  m_background   = m_config->get_number(prefix + "background", "kelvin");
  m_excess       = m_config->get_number(prefix + "excess", "kelvin");
  m_scale_height = m_config->get_number(prefix + "scale_height", "m");
  m_cos          = m_config->get_number(prefix + "cos", "kelvin");
  m_cos_gradient = m_config->get_number(prefix + "cos_gradient", "kelvin m^-1");
  m_sin          = m_config->get_number(prefix + "sin", "kelvin");
  m_sin_gradient = m_config->get_number(prefix + "sin_gradient", "kelvin m^-1");
  m_z_min        = m_config->get_number(prefix + "z_min", "m");
  m_z_max        = m_config->get_number(prefix + "z_max", "m");

  if (not(m_scale_height > 0.0)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "%sscale_height has to be positive (got %f m)",
                                  prefix.c_str(), m_scale_height);
  }

  if (not(m_z_min < m_z_max)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "%sz_min (%f m) has to be below %sz_max (%f m)",
                                  prefix.c_str(), m_z_min, prefix.c_str(), m_z_max);
  }

  m_temperature = allocate_temperature(grid);

  m_surface.metadata(0)
      .long_name("ice upper surface elevation")
      .units("m")
      .standard_name("surface_altitude");
}

void CARRA2Bias::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing the elevation- and season-dependent air temperature bias "
                 "correction 'carra2_bias'...\n"
                 "  T = T_input + b(z, t), with z the surface elevation and t the fraction of the year:\n"
                 "  b = %.2f K + %.2f K exp(-z / %.0f m)\n"
                 "      + (%.2f K + %.3f K/km z) cos(2 pi t) + (%.2f K + %.3f K/km z) sin(2 pi t),\n"
                 "  with z held at %.0f m below and at %.0f m above these elevations.\n",
                 m_background, m_excess, m_scale_height,
                 m_cos, m_cos_gradient * 1000.0, m_sin, m_sin_gradient * 1000.0,
                 m_z_min, m_z_max);
}

double CARRA2Bias::bias(double z, double cos_t, double sin_t) const {
  // hold the elevation at the thresholds outside [z_min, z_max]
  double Z = std::min(std::max(z, m_z_min), m_z_max);

  return (m_background + m_excess * std::exp(-Z / m_scale_height) +
          (m_cos + m_cos_gradient * Z) * cos_t +
          (m_sin + m_sin_gradient * Z) * sin_t);
}

void CARRA2Bias::mean_cycle(double t, double dt, double &cos_t, double &sin_t) const {
  // Sample the step at (at least) daily resolution: a step can span a season or more,
  // and the seasonal terms average to zero over a whole year.
  const double day = 86400.0;
  const int N      = std::max(1, std::min(10000, static_cast<int>(std::ceil(dt / day))));

  cos_t = 0.0;
  sin_t = 0.0;
  for (int k = 0; k < N; ++k) {
    double phase = 2.0 * M_PI * time().year_fraction(t + (k + 0.5) * dt / N);
    cos_t += std::cos(phase);
    sin_t += std::sin(phase);
  }
  cos_t /= N;
  sin_t /= N;
}

void CARRA2Bias::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  // keep a copy of the surface elevation for temp_time_series_impl()
  m_surface.copy_from(geometry.ice_surface_elevation);

  double cos_t = 0.0, sin_t = 0.0;
  mean_cycle(t, dt, cos_t, sin_t);

  m_temperature->copy_from(m_input_model->air_temperature());

  auto &T = *m_temperature;

  array::AccessScope list{ &T, &m_surface };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    T(i, j) += bias(m_surface(i, j), cos_t, sin_t);
  }
}

const array::Scalar &CARRA2Bias::air_temperature_impl() const {
  return *m_temperature;
}

void CARRA2Bias::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();
  m_surface.begin_access();
}

void CARRA2Bias::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();
  m_surface.end_access();
}

void CARRA2Bias::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_cos_values.resize(ts.size());
  m_sin_values.resize(ts.size());

  for (size_t k = 0; k < ts.size(); ++k) {
    double phase    = 2.0 * M_PI * time().year_fraction(ts[k]);
    m_cos_values[k] = std::cos(phase);
    m_sin_values[k] = std::sin(phase);
  }
}

void CARRA2Bias::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->temp_time_series(i, j, result);

  double z = m_surface(i, j);

  for (size_t k = 0; k < result.size(); ++k) {
    result[k] += bias(z, m_cos_values[k], m_sin_values[k]);
  }
}

} // end of namespace atmosphere
} // end of namespace pism
