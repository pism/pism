// Copyright (C) 2004-2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "IceModel.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/pism_signal.h"

namespace pism {

//! Catch signals -USR1, -USR2 and -TERM.
/*!
Signal `SIGTERM` makes PISM end, saving state under original `-o` name
(or default name).  We also add an indication to the history attribute
of the output NetCDF file.

Signal `SIGUSR1` makes PISM save state under a filename based on the
the name of the executable (e.g. `pismr` or `pismv`) and the current
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
                                 m_time->date().c_str()));
    // Tell the caller that the user requested an early termination of
    // the run.
    return 1;
  }

  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "pism-%s.nc",
             m_time->date().c_str());
    m_log->message(1,
       "\ncaught signal SIGUSR1:  Writing intermediate file `%s' and flushing time series.\n\n",
       file_name);
    pism_signal = 0;

    PIO file(m_grid->com, m_config->get_string("output.format"), file_name, PISM_READWRITE_MOVE);
    save_variables(file, INCLUDE_MODEL_STATE, m_output_vars);

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


void IceModel::update_run_stats() {

  // timing stats
  // MYPPH stands for "model years per processor hour"
  double
    wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time),
    proc_hours       = m_grid->size() * wall_clock_hours,
    mypph            = units::convert(m_sys,
                                      m_time->current() - m_time->start(),
                                      "seconds", "years") / proc_hours;

  // time-independent info
  {
    m_run_stats.set_string("source", std::string("PISM ") + PISM_Revision);
    m_run_stats.set_string("long_name", "Run statistics");
  }

  m_run_stats.set_double("wall_clock_hours", wall_clock_hours);
  m_run_stats.set_double("processor_hours", proc_hours);
  m_run_stats.set_double("model_years_per_processor_hour", mypph);
}

//! Get time and user/host name and add it to the given string.
void  IceModel::prepend_history(const std::string &str) {
  m_output_global_attributes.set_string("history",
                                        username_prefix(m_grid->com) + (str + "\n") +
                                        m_output_global_attributes.get_string("history"));
}

//! Check if the thickness of the ice is too large.
/*! Return true if the ice thickness exceeds the height of the computational domain.
 */
bool check_maximum_ice_thickness(const IceModelVec2S &ice_thickness) {
  IceGrid::ConstPtr grid = ice_thickness.grid();

  const double Lz = grid->Lz();

  IceModelVec::AccessList list(ice_thickness);

  unsigned int counter = 0;
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_thickness(i, j) > Lz) {
      counter += 1;
    }
  }

  if (GlobalSum(grid->com, counter) > 0) {
    return true;
  } else {
    return false;
  }
}

//! Return the grid used by this model.
IceGrid::Ptr IceModel::grid() const {
  return m_grid;
}

//! Return the context this model is running in.
Context::Ptr IceModel::ctx() const {
  return m_ctx;
}

} // end of namespace pism
