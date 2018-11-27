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

#include "pism/util/pism_utilities.hh"
#include "pism/util/Profiling.hh"

namespace pism {

//! Computes the maximum time-step we can take and still hit all `-save_times`.
MaxTimestep IceModel::save_max_timestep(double my_t) {

  if ((not m_save_snapshots) or
      (not m_config->get_boolean("time_stepping.hit_save_times"))) {
    return MaxTimestep("reporting (-save_times)");
  }

  return reporting_max_timestep(m_snapshot_times, my_t, "reporting (-save_times)");
}

//! Initializes the snapshot-saving mechanism.
void IceModel::init_snapshots() {
  m_current_snapshot = 0;

  m_snapshots_filename = m_config->get_string("output.snapshot.file");
  bool filename_set = not m_snapshots_filename.empty();

  auto save_times = m_config->get_string("output.snapshot.times");
  bool times_set = not save_times.empty();

  bool split = m_config->get_boolean("output.snapshot.split");

  m_snapshot_vars = output_variables(m_config->get_string("output.snapshot.size"));

  if (filename_set ^ times_set) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "you need to set both output.snapshot.file and output.snapshot.times"
                       " to save snapshots.");
  }

  if (not filename_set and not times_set) {
    m_save_snapshots = false;
    return;
  }

  try {
    // parse
    std::vector<double> times = m_time->parse_times(save_times);

    // discard times before the beginning and after the end of the run
    m_snapshot_times.clear();
    for (const auto &t : times) {
      if (t >= m_time->start() and t <= m_time->end()) {
        m_snapshot_times.push_back(t);
      }
    }
  } catch (RuntimeError &e) {
    e.add_context("processing output.snapshot.times");
    throw;
  }

  if (m_snapshot_times.size() == 0) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "output.snapshot.times was set, but all requested times"
                       " are outside of the modeled time interval");
  }

  m_save_snapshots = true;
  m_snapshots_file_is_ready = false;
  m_split_snapshots = false;

  if (split) {
    m_split_snapshots = true;
  } else if (not ends_with(m_snapshots_filename, ".nc")) {
    m_log->message(2,
               "PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
  }

  if (split) {
    m_log->message(2, "saving snapshots to '%s+year.nc'; ",
               m_snapshots_filename.c_str());
  } else {
    m_log->message(2, "saving snapshots to '%s'; ",
               m_snapshots_filename.c_str());
  }

  m_log->message(2, "times requested: %s\n", save_times.c_str());
}

  //! Writes a snapshot of the model state (if necessary)
void IceModel::write_snapshot() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
  // value is never used, because saving_after
  // is only used if save_now == true, and in
  // this case saving_after is guaranteed to be
  // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];

  // determine if the user set the -save_times and -save_file options
  if (not m_save_snapshots) {
    return;
  }

  // do we need to save *now*?
  if ((m_time->current() >= m_snapshot_times[m_current_snapshot]) and
      (m_current_snapshot < m_snapshot_times.size())) {
    saving_after = m_snapshot_times[m_current_snapshot];

    while ((m_current_snapshot < m_snapshot_times.size()) and
           (m_snapshot_times[m_current_snapshot] <= m_time->current())) {
      m_current_snapshot++;
    }
  } else {
    // we don't need to save now, so just return
    return;
  }

  // flush time-series buffers
  flush_timeseries();

  if (m_split_snapshots) {
    m_snapshots_file_is_ready = false;    // each snapshot is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s_%s.nc",
             m_snapshots_filename.c_str(), m_time->date(saving_after).c_str());
  } else {
    strncpy(filename, m_snapshots_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  m_log->message(2,
             "saving snapshot to %s at %s, for time-step goal %s\n",
             filename, m_time->date().c_str(),
             m_time->date(saving_after).c_str());

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("io.snapshots");
  IO_Mode mode = m_snapshots_file_is_ready ? PISM_READWRITE : PISM_READWRITE_MOVE;
  {
    PIO file(m_grid->com, m_config->get_string("output.format"), filename, mode);

    if (not m_snapshots_file_is_ready) {
      write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

      m_snapshots_file_is_ready = true;
    }

    write_run_stats(file);

    save_variables(file, INCLUDE_MODEL_STATE, m_snapshot_vars);
  }
  profiling.end("io.snapshots");
}

} // end of namespace pism
