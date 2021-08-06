/* Copyright (C) 2018, 2019, 2020, 2021 PISM Authors
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

#include "ScalarForcing.hh"

#include <cassert>              // assert()
#include <cmath>                // std::floor()

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Timeseries.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/io/File.hh"

namespace pism {

ScalarForcing::ScalarForcing(const Context &ctx,
                             const std::string &prefix,
                             const std::string &variable_name,
                             const std::string &units,
                             const std::string &glaciological_units,
                             const std::string &long_name)
  : m_period(0.0), m_period_start(0.0), m_current(0.0) {

  Config::ConstPtr config = ctx.config();

  m_data.reset(new Timeseries(ctx.unit_system(), variable_name));
  m_data->variable().set_string("units", units);
  m_data->variable().set_string("glaciological_units", glaciological_units);
  m_data->variable().set_string("long_name", long_name);

  // Read data from a NetCDF file
  {
    auto filename = config->get_string(prefix + ".file");

    if (filename.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "%s.file is required", prefix.c_str());
    }

    ctx.log()->message(2,
                       "  reading %s data from forcing file %s...\n",
                       m_data->name().c_str(), filename.c_str());

    try {
      File file(ctx.com(), filename, PISM_NETCDF3, PISM_READONLY);
      {
        auto time_units = ctx.time()->units_string();
        m_data->read(file, time_units, *ctx.log());
      }
    } catch (RuntimeError &e) {
      e.add_context("while reading %s (%s) from '%s'",
                    long_name.c_str(), variable_name.c_str(), filename.c_str());
      throw;
    }

    bool periodic = config->get_flag(prefix + ".periodic");

    if (periodic) {
      auto T = m_data->time_interval();

      m_period = T[1] - T[0];
      assert(m_period > 0.0);

      m_period_start = T[0];
    } else {
      m_period = 0.0;
      m_period_start = 0.0;
    }

  }
}

void ScalarForcing::update(double t, double dt) {
  m_current = value(t + 0.5 * dt);
}

double ScalarForcing::value() const {
  return m_current;
}

double ScalarForcing::value(double t) const {
  if (m_period > 0.0) {
    double T = t - m_period_start;

    T -= std::floor(T / m_period) * m_period;

    return (*m_data)(m_period_start + T);
  }

  return (*m_data)(t);
}

} // end of namespace pism
