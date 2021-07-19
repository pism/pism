/* Copyright (C) 2018, 2019, 2020 PISM Authors
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

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Timeseries.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/io/File.hh"

namespace pism {

ScalarForcing::ScalarForcing(std::shared_ptr<const Context> ctx,
                             const std::string &prefix,
                             const std::string &variable_name,
                             const std::string &units,
                             const std::string &glaciological_units,
                             const std::string &long_name)
  : m_ctx(ctx), m_period(0), m_reference_time(0.0), m_current(0.0) {

  Config::ConstPtr config = ctx->config();

  m_prefix = prefix;

  m_data.reset(new Timeseries(ctx->com(), ctx->unit_system(), variable_name));
  m_data->variable().set_string("units", units);
  m_data->variable().set_string("glaciological_units", glaciological_units);
  m_data->variable().set_string("long_name", long_name);
}

ScalarForcing::~ScalarForcing() {
  // empty
}

void ScalarForcing::init() {

  Config::ConstPtr config = m_ctx->config();

  auto   filename       = config->get_string(m_prefix + ".file");
  int    period         = config->get_number(m_prefix + ".period");
  double reference_year = config->get_number(m_prefix + ".reference_year");

  if (filename.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "%s.file is required", m_prefix.c_str());
  }

  if (period < 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid %s.period %d (period length cannot be negative)",
                                  m_prefix.c_str(), period);
  }

  m_period = (unsigned int)period;

  m_reference_time = units::convert(m_ctx->unit_system(),
                                    reference_year, "years", "seconds");

  m_ctx->log()->message(2,
                        "  reading %s data from forcing file %s...\n",
                        m_data->name().c_str(), filename.c_str());

  File file(m_ctx->com(), filename, PISM_NETCDF3, PISM_READONLY);
  {
    m_data->read(file, *m_ctx->time(), *m_ctx->log());
  }
}

void ScalarForcing::update(double t, double dt) {
  m_current = value(t + 0.5 * dt);
}

double ScalarForcing::value() const {
  return m_current;
}

double ScalarForcing::value(double t) const {
  t = m_ctx->time()->modulo(t - m_reference_time, m_period);

  return (*m_data)(t);
}

} // end of namespace pism
