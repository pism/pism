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

#include "ScalarForcing.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Timeseries.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/io/PIO.hh"

namespace pism {

ScalarForcing::ScalarForcing(Context::ConstPtr ctx,
                             const std::string &option_prefix,
                             const std::string &offset_name,
                             const std::string &units,
                             const std::string &glaciological_units,
                             const std::string &long_name)
  : m_ctx(ctx), m_bc_period(0), m_bc_reference_time(0.0), m_current(0.0) {

  Config::ConstPtr config = ctx->config();

  m_option_prefix = option_prefix;

  m_offset.reset(new Timeseries(ctx->com(), ctx->unit_system(),
                                offset_name,
                                config->get_string("time.dimension_name")));
  m_offset->variable().set_string("units", units);
  m_offset->variable().set_string("glaciological_units", glaciological_units);
  m_offset->variable().set_string("long_name", long_name);

  m_offset->dimension().set_string("units", ctx->time()->units_string());
}

ScalarForcing::~ScalarForcing() {
  // empty
}

void ScalarForcing::init() {
  options::String file(m_option_prefix + "_file", "Specifies a file with scalar offsets");
  options::Integer period(m_option_prefix + "_period",
                          "Specifies the length of the climate data period", 0);
  options::Real bc_reference_year(m_option_prefix + "_reference_year",
                                  "Boundary condition reference year", 0.0);

  if (not file.is_set()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "command-line option %s_file is required.",
                                  m_option_prefix.c_str());
  }

  if (period.value() < 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid %s_period %d (period length cannot be negative)",
                                  m_option_prefix.c_str(), period.value());
  }

  m_bc_period = (unsigned int)period;

  if (bc_reference_year.is_set()) {
    m_bc_reference_time = units::convert(m_ctx->unit_system(),
                                         bc_reference_year, "years", "seconds");
  } else {
    m_bc_reference_time = 0;
  }

  m_ctx->log()->message(2,
                        "  reading %s data from forcing file %s...\n",
                        m_offset->name().c_str(), file->c_str());

  PIO nc(m_ctx->com(), "netcdf3", file, PISM_READONLY);
  {
    m_offset->read(nc, *m_ctx->time(), *m_ctx->log());
  }
}

void ScalarForcing::update(double t, double dt) {
  t = m_ctx->time()->mod(t - m_bc_reference_time, m_bc_period);

  m_current = (*m_offset)(t + 0.5 * dt);
}

double ScalarForcing::value() const {
  return m_current;
}

} // end of namespace pism
