// Copyright (C) 2004-2019, 2021, 2023, 2024 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include <petscsys.h>

#include "pism/icemodel/IceModel.hh"

#include "pism/util/Grid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_signal.h"

namespace pism {

//! Catch signals -USR1, -USR2 and -TERM.
/*!
Signal `SIGTERM` makes PISM end, saving state under original `-o` name
(or default name).  We also add an indication to the history attribute
of the output NetCDF file.

Signal `SIGUSR1` makes PISM save state under a filename based on the
the name of the executable (e.g. `pism`) and the current
model year.  In addition the time series (`-ts_file`, etc.) is flushed out
There is no indication of these actions in the history attribute of the output (`-o`)
NetCDF file because there is no effect on it, but there is an indication at `stdout`.

Signal `SIGUSR2` makes PISM flush time-series, without saving model state.
 */
int IceModel::process_signals() {

  if (pism_signal == SIGTERM) {
    m_log->message(1,
       "\ncaught signal SIGTERM:  EXITING EARLY and saving with original filename.\n");

    prepend_history(pism::printf("EARLY EXIT caused by signal SIGTERM. Completed timestep at time=%s.",
                                 m_time->date(m_time->current()).c_str()));
    // Tell the caller that the user requested an early termination of
    // the run.
    return 1;
  }

  if (pism_signal == SIGUSR1) {
    auto date_without_spaces = replace_character(m_time->date(m_time->current()), ' ', '_');
    auto file_name = pism::printf("pism-%s.nc", date_without_spaces.c_str());
    m_log->message(1,
       "\ncaught signal SIGUSR1:  Writing intermediate file `%s' and flushing time series.\n\n",
                   file_name.c_str());
    pism_signal = 0;

    File file(m_grid->com,
              file_name,
              string_to_backend(m_config->get_string("output.format")),
              io::PISM_READWRITE_MOVE);
    save_variables(file, INCLUDE_MODEL_STATE, m_output_vars, m_time->current());

    // flush all the time-series buffers:
    flush_timeseries();
  }

  if (pism_signal == SIGUSR2) {
    m_log->message(1,
       "\ncaught signal SIGUSR2:  Flushing time series.\n\n");
    pism_signal = 0;

    // flush all the time-series buffers:
    flush_timeseries();
  }

  return 0;
}


VariableMetadata IceModel::run_stats() const {

  VariableMetadata result("run_stats", m_sys);

  result["source"]    = std::string("PISM ") + pism::revision;
  result["long_name"] = "Run statistics";

  // timing stats
  double
    wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time),
    proc_hours       = m_grid->size() * wall_clock_hours,
    model_years      = units::convert(m_sys, m_time->current() - m_time->start(),
                                      "seconds", "years");

  result["wall_clock_hours"]               = { wall_clock_hours };
  result["processor_hours"]                = { proc_hours };
  result["model_years_per_processor_hour"] = { model_years / proc_hours };
  result["number_of_time_steps"]           = { (double)m_step_counter };

  return result;
}

//! Get time and user/host name and add it to the given string.
void  IceModel::prepend_history(const std::string &str) {
  m_output_global_attributes.set_string("history",
                                        username_prefix(m_grid->com) + (str + "\n") +
                                        m_output_global_attributes.get_string("history"));
}

//! Return the grid used by this model.
std::shared_ptr<Grid> IceModel::grid() const {
  return m_grid;
}

//! Return the context this model is running in.
std::shared_ptr<Context> IceModel::ctx() const {
  return m_ctx;
}

} // end of namespace pism
