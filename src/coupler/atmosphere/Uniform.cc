/* Copyright (C) 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "Uniform.hh"

#include "pism/geometry/Geometry.hh"

namespace pism {
namespace atmosphere {

Uniform::Uniform(IceGrid::ConstPtr grid)
  : AtmosphereModel(grid, nullptr) {
  m_precipitation = allocate_precipitation(grid);
  m_temperature = allocate_temperature(grid);
}

void Uniform::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the test atmosphere model...\n");

  m_temperature->set(m_config->get_double("atmosphere.uniform.temperature", "Kelvin"));
  m_precipitation->set(m_config->get_double("atmosphere.uniform.precipitation", "kg m-2 s-1"));
}

void Uniform::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) t;
  (void) dt;
}

const IceModelVec2S& Uniform::mean_precipitation_impl() const {
  return *m_precipitation;
}

const IceModelVec2S& Uniform::mean_annual_temp_impl() const {
  return *m_temperature;
}

void Uniform::begin_pointwise_access_impl() const {
  m_precipitation->begin_access();
  m_temperature->begin_access();
}

void Uniform::end_pointwise_access_impl() const {
  m_precipitation->end_access();
  m_temperature->end_access();
}


void Uniform::temp_time_series_impl(int i, int j, std::vector<double> &values) const {
  for (size_t k = 0; k < m_ts_times.size(); ++k) {
    values[k] = (*m_temperature)(i, j);
  }
}

void Uniform::precip_time_series_impl(int i, int j, std::vector<double> &values) const {
  for (size_t k = 0; k < m_ts_times.size(); ++k) {
    values[k] = (*m_precipitation)(i, j);
  }
}

} // end of namespace atmosphere
} // end of namespace pism
