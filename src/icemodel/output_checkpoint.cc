/* Copyright (C) 2017, 2019, 2022, 2023, 2024, 2025 PISM Authors
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
    m_checkpoint_filename = filename_add_suffix(m_output_filename, "_checkpoint", "");
  }

  m_checkpoint_vars = output_variables(m_config->get_string("output.checkpoint.size"));
  m_last_checkpoint_time = 0.0;

  {
    m_checkpoint_file_contents = pism::combine(common_metadata(), state_variables());
    m_checkpoint_file_contents =
        pism::combine(m_checkpoint_file_contents, diagnostic_variables(m_checkpoint_vars));
  }
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
    OutputFile file(m_output_writer, m_checkpoint_filename);
    // Ensure that the checkpoint file is closed to force PISM to open a new file every
    // time we write a checkpoint, moving the old file aside if it exists.
    file.close();

    {
      // define time dimension *without* time bounds
      define_time(file);
      define_variables(file, m_checkpoint_file_contents);
    }

    {
      write_config(*m_config, "pism_config", file);
      file.append_time(m_time->current());
      write_state(file);
      write_diagnostics(file, m_checkpoint_vars);
      write_run_stats(file);
    }
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
