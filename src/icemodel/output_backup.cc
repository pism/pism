/* Copyright (C) 2017 PISM Authors
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

//! Initialize the backup (snapshot-on-wallclock-time) mechanism.
void IceModel::init_backups() {

  std::string backup_file = m_config->get_string("output.file_name");
  if (not backup_file.empty()) {
    m_backup_filename = filename_add_suffix(backup_file, "_backup", "");
  } else {
    m_backup_filename = "pism_backup.nc";
  }

  m_backup_vars = output_variables(m_config->get_string("output.backup_size"));
  m_last_backup_time = 0.0;
}

  //! Write a backup (i.e. an intermediate result of a run).
void IceModel::write_backup() {

  double backup_interval = m_config->get_double("output.backup_interval");

  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time);

  if (wall_clock_hours - m_last_backup_time < backup_interval) {
    return;
  }

  const Profiling &profiling = m_ctx->profiling();

  m_last_backup_time = wall_clock_hours;

  // create a history string:

  m_log->message(2,
                 "  [%s] Saving an automatic backup to '%s' (%1.3f hours after the beginning of the run)\n",
                 timestamp(m_grid->com).c_str(), m_backup_filename.c_str(), wall_clock_hours);

  double backup_start_time = get_time();
  profiling.begin("io.backup");
  {
    PIO file(m_grid->com, m_config->get_string("output.format"),
             m_backup_filename, PISM_READWRITE_MOVE);

    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);
    write_run_stats(file);

    save_variables(file, INCLUDE_MODEL_STATE, m_backup_vars);
  }
  profiling.end("io.backup");
  double backup_end_time = get_time();

  // Also flush time-series:
  flush_timeseries();

  m_log->message(2,
                 "  [%s] Done saving an automatic backup in %f seconds (%f minutes).\n",
                 timestamp(m_grid->com).c_str(),
                 backup_end_time - backup_start_time,
                 (backup_end_time - backup_start_time) / 60.0);

}

} // end of namespace pism
