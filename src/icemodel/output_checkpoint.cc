/* Copyright (C) 2017, 2019, 2022, 2023 PISM Authors
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

namespace pism {

//! Initialize checkpointing (snapshot-on-wallclock-time) mechanism.
void IceModel::init_checkpoints() {

  m_checkpoint_filename = m_config->get_string("output.checkpoint.file");

  if (m_checkpoint_filename.empty()) {
    std::string output_file = m_config->get_string("output.file");
    if (not output_file.empty()) {
      m_checkpoint_filename = filename_add_suffix(output_file, "_checkpoint", "");
    } else {
      m_checkpoint_filename = "pism_checkpoint.nc";
    }
  }

  m_checkpoint_vars = output_variables(m_config->get_string("output.checkpoint.size"));
  m_last_checkpoint_time = 0.0;
}

//! Write a checkpoint (i.e. an intermediate result of a run).
/*!
 * Returns `true` if PISM has to stop, `false` otherwise.
 */
bool IceModel::write_checkpoint() {

  double checkpoint_interval = m_config->get_number("output.checkpoint.interval");

  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time);

  if (wall_clock_hours - m_last_checkpoint_time < checkpoint_interval) {
    return false;
  }

  const Profiling &profiling = m_ctx->profiling();

  m_last_checkpoint_time = wall_clock_hours;

  // create a history string:

  m_log->message(2,
                 "  [%s] Saving a checkpoint to '%s' (%1.3f hours after the beginning of the run)\n",
                 timestamp(m_grid->com).c_str(), m_checkpoint_filename.c_str(), wall_clock_hours);

  double checkpoint_start_time = get_time(m_grid->com);
  profiling.begin("io.checkpoint");
  {
    File file(m_grid->com,
              m_checkpoint_filename,
              string_to_backend(m_config->get_string("output.format")),
              io::PISM_READWRITE_MOVE,
              m_ctx->pio_iosys_id());

    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);
    write_run_stats(file, run_stats());

    save_variables(file, INCLUDE_MODEL_STATE, m_checkpoint_vars, m_time->current());
  }
  profiling.end("io.checkpoint");
  double checkpoint_end_time = get_time(m_grid->com);

  // Also flush time-series:
  flush_timeseries();

  m_log->message(2,
                 "  [%s] Done saving a checkpoint in %f seconds (%f minutes).\n",
                 timestamp(m_grid->com).c_str(),
                 checkpoint_end_time - checkpoint_start_time,
                 (checkpoint_end_time - checkpoint_start_time) / 60.0);

  return m_config->get_flag("output.checkpoint.exit");
}

} // end of namespace pism
