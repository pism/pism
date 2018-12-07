/* Copyright (C) 2017, 2018 PISM Authors
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

#include "IceModel.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

//! Initializes the code writing scalar time-series.
void IceModel::init_timeseries() {

  m_ts_filename = m_config->get_string("output.timeseries.filename");

  auto times = m_config->get_string("output.timeseries.times");
  bool times_set = not times.empty();

  if (times_set xor not m_ts_filename.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "you need to specity both -ts_file and -ts_times"
                       " to save scalar diagnostic time-series.");
  }

  if (m_ts_filename.empty()) {
    return;
  }

  try {
    *m_ts_times = m_time->parse_times(times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -ts_times argument %s", times.c_str());
    throw;
  }

  m_log->message(2, "  saving scalar time-series to '%s'\n", m_ts_filename.c_str());
  m_log->message(2, "  times requested: %s\n", times.c_str());

  auto vars = set_split(m_config->get_string("output.timeseries.variables"), ',');
  if (not vars.empty()) {
    m_log->message(2, "variables requested: %s\n", set_join(vars, ",").c_str());
  }

  // prepare the output file
  {
    // default behavior is to move the file aside if it exists already; option allows appending
    bool append = m_config->get_boolean("output.timeseries.append");
    IO_Mode mode = append ? PISM_READWRITE : PISM_READWRITE_MOVE;
    PIO file(m_grid->com, "netcdf3", m_ts_filename, mode);      // Use NetCDF-3 to write time-series.
    // add the last saved time to the list of requested times so that the first time is interpreted
    // as the end of a reporting time step
    if (append and file.inq_dimlen("time") > 0) {
      double
        epsilon = 1.0,          // one second
        t       = 0.0;
      file.inq_dim_limits("time", NULL, &t);
      // add this time only if it is strictly before the first requested one
      if (t + epsilon < m_ts_times->front()) {
        m_ts_times->insert(m_ts_times->begin(), t);
      }
    }

    write_metadata(file, SKIP_MAPPING, PREPEND_HISTORY);
    write_run_stats(file);

    // initialize scalar diagnostics
    for (auto d : m_ts_diagnostics) {
      d.second->init(file, m_ts_times);
    }
  }
}

//! Computes the maximum time-step we can take and still hit all `-ts_times`.
MaxTimestep IceModel::ts_max_timestep(double my_t) {

  if ((not m_config->get_boolean("time_stepping.hit_ts_times")) or
      m_ts_diagnostics.empty()) {
    return MaxTimestep("reporting (-ts_times)");
  }

  return reporting_max_timestep(*m_ts_times, my_t, "reporting (-ts_times)");
}

//! Flush scalar time-series.
void IceModel::flush_timeseries() {
  // flush all the time-series buffers:
  for (auto d : m_ts_diagnostics) {
    d.second->flush();
  }

  // update run_stats in the time series output file
  if (not m_ts_diagnostics.empty()) {
    PIO file(m_grid->com, "netcdf3", m_ts_filename, PISM_READWRITE);
    write_run_stats(file);
  }
}

} // end of namespace pism
