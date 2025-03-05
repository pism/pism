/* Copyright (C) 2017, 2018, 2019, 2021, 2023, 2024 PISM Authors
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

#include "pism/icemodel/IceModel.hh"

#include "pism/util/pism_utilities.hh"
#include "pism/util/Profiling.hh"
#include <memory>

namespace pism {

//! Computes the maximum time-step we can take and still hit all `-save_times`.
MaxTimestep IceModel::save_max_timestep(double my_t) {

  if (m_snapshots_filename.empty() or (not m_config->get_flag("time_stepping.hit_save_times"))) {
    return MaxTimestep("reporting (-save_times)");
  }

  double eps = m_config->get_number("time_stepping.resolution");

  return reporting_max_timestep(m_snapshot_times, my_t, eps,
                                "reporting (-save_times)");
}

//! Initializes the snapshot-saving mechanism.
void IceModel::init_snapshots() {
  m_current_snapshot = 0;

  m_snapshots_filename = m_config->get_string("output.snapshot.file");
  auto save_times      = m_config->get_string("output.snapshot.times");
  m_snapshot_vars      = output_variables(m_config->get_string("output.snapshot.size"));
  m_split_snapshots    = m_config->get_flag("output.snapshot.split");

  {
    bool filename_set = not m_snapshots_filename.empty();
    bool times_set    = not save_times.empty();

    if (filename_set ^ times_set) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "you need to set both output.snapshot.file and output.snapshot.times"
                         " to save snapshots.");
    }

    if (not (filename_set and times_set)) {
      return;
    }
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

  if (m_snapshot_times.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "output.snapshot.times was set, but all requested times"
                       " are outside of the modeled time interval");
  }

  if (m_split_snapshots) {
    m_log->message(2, "saving snapshots to '%s+year.nc'; ", m_snapshots_filename.c_str());
  } else {
    m_log->message(2, "saving snapshots to '%s'; ", m_snapshots_filename.c_str());

    if (not ends_with(m_snapshots_filename, ".nc")) {
      m_log->message(2, "PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
    }
  }

  m_log->message(2, "times requested: %s\n", save_times.c_str());
}

//! Writes a snapshot of the model state (if necessary)
void IceModel::write_snapshot() {
  // initialize to avoid compiler warning; this value is never used, because saving_after
  // is only used if save_now == true, and in this case saving_after is guaranteed to be
  // initialized. See the code below.
  double saving_after = -1.0e30;

  // determine if the user set the -save_times and -save_file options
  if (m_snapshots_filename.empty()) {
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

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("io.snapshots");
  std::string filename;
  if (m_snapshot_file == nullptr) {
    if (m_split_snapshots) {
      auto date_without_spaces  = replace_character(m_time->date(saving_after), ' ', '_');
      filename =
          pism::printf("%s_%s.nc", m_snapshots_filename.c_str(), date_without_spaces.c_str());
    } else {
      filename = m_snapshots_filename;
    }

    m_snapshot_file = std::make_shared<File>(
        m_grid->com, filename, string_to_backend(m_config->get_string("output.format")),
        io::PISM_READWRITE_MOVE);

    write_metadata(*m_snapshot_file, WRITE_MAPPING, PREPEND_HISTORY);
  }

  {
    m_log->message(2, "saving snapshot to %s at %s, for time-step goal %s\n", filename.c_str(),
                   m_time->date(m_time->current()).c_str(), m_time->date(saving_after).c_str());
    write_run_stats(*m_snapshot_file, run_stats());
    save_variables(*m_snapshot_file, INCLUDE_MODEL_STATE, m_snapshot_vars, m_time->current());
  }

  if (m_split_snapshots) {
    m_snapshot_file.reset();
  }

  profiling.end("io.snapshots");
}

} // end of namespace pism
