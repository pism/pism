// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <sstream>
#include <cstring>
#include <petscsys.h>

#include "iceModel.hh"

#include "base/energy/bedrockThermalUnit.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "coupler/PISMSurface.hh"
#include "enthalpyConverter.hh"
#include "pism_signal.h"
#include "base/util/pism_utilities.hh"

namespace pism {


//! Virtual.  Does nothing in `IceModel`.  Derived classes can do more computation in each time step.
void IceModel::additionalAtStartTimestep() {
  // empty
}


//! Virtual.  Does nothing in `IceModel`.  Derived classes can do more computation in each time step.
void IceModel::additionalAtEndTimestep() {
  // empty
}

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
int IceModel::endOfTimeStepHook() {

  if (pism_signal == SIGTERM) {
    m_log->message(1,
       "\ncaught signal SIGTERM:  EXITING EARLY and saving with original filename.\n");
    char str[TEMPORARY_STRING_LENGTH];
    snprintf(str, sizeof(str),
       "EARLY EXIT caused by signal SIGTERM.  Completed timestep at time=%s.",
             m_time->date().c_str());
    stampHistory(str);
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
    dumpToFile(file_name);

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


//! Build a history string from the command which invoked PISM.
void  IceModel::stampHistoryCommand() {

  char startstr[TEMPORARY_STRING_LENGTH];

  snprintf(startstr, sizeof(startstr),
           "PISM (%s) started on %d procs.", PISM_Revision, (int)m_grid->size());
  stampHistory(std::string(startstr));

  m_output_global_attributes.set_string("history",
                               pism_args_string() + m_output_global_attributes.get_string("history"));
}

void IceModel::update_run_stats() {
  PetscErrorCode ierr;

  // timing stats
  // MYPPH stands for "model years per processor hour"
  double
    wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time),
    proc_hours       = m_grid->size() * wall_clock_hours,
    mypph            = units::convert(m_sys,
                                      m_time->current() - m_time->start(),
                                      "seconds", "years") / proc_hours;

  // get PETSc's reported number of floating point ops (*not* per time) on this
  //   process, then sum over all processes
  PetscLogDouble my_flops = 0.0;
  ierr = PetscGetFlops(&my_flops);
  PISM_CHK(ierr, "PetscGetFlops");

  double flops = GlobalSum(m_grid->com, my_flops);

  run_stats.set_double("wall_clock_hours", wall_clock_hours);
  run_stats.set_double("processor_hours", proc_hours);
  run_stats.set_double("model_years_per_processor_hour", mypph);
  run_stats.set_double("PETSc_MFlops", flops * 1.0e-6);
  run_stats.set_double("grid_dx_meters", m_grid->dx());
  run_stats.set_double("grid_dy_meters", m_grid->dy());
  run_stats.set_double("grid_dz_min_meters", m_grid->dz_min());
  run_stats.set_double("grid_dz_max_meters", m_grid->dz_max());
  if (btu != NULL) {
    run_stats.set_double("grid_dzb_meters", btu->vertical_spacing());
  }
  run_stats.set_string("source", std::string("PISM ") + PISM_Revision);

  run_stats.set_double("grounded_basal_ice_flux_cumulative", grounded_basal_ice_flux_cumulative);
  run_stats.set_double("nonneg_rule_flux_cumulative", nonneg_rule_flux_cumulative);
  run_stats.set_double("sub_shelf_ice_flux_cumulative", sub_shelf_ice_flux_cumulative);
  run_stats.set_double("surface_ice_flux_cumulative", surface_ice_flux_cumulative);
  run_stats.set_double("sum_divQ_SIA_cumulative", sum_divQ_SIA_cumulative);
  run_stats.set_double("sum_divQ_SSA_cumulative", sum_divQ_SSA_cumulative);
  run_stats.set_double("Href_to_H_flux_cumulative", Href_to_H_flux_cumulative);
  run_stats.set_double("H_to_Href_flux_cumulative", H_to_Href_flux_cumulative);
  run_stats.set_double("discharge_flux_cumulative", discharge_flux_cumulative);
}

//! Build the particular history string associated to the end of a PISM run,
//! including a minimal performance assessment.
void  IceModel::stampHistoryEnd() {

  update_run_stats();

  // build and put string into global attribute "history"
  char str[TEMPORARY_STRING_LENGTH];

  snprintf(str, TEMPORARY_STRING_LENGTH,
    "PISM done.  Performance stats: %.4f wall clock hours, %.4f proc.-hours, %.4f model years per proc.-hour, PETSc MFlops = %.2f.",
           run_stats.get_double("wall_clock_hours"),
           run_stats.get_double("processor_hours"),
           run_stats.get_double("model_years_per_processor_hour"),
           run_stats.get_double("PETSc_MFlops"));

  stampHistory(str);
}


//! Get time and user/host name and add it to the given string.
void  IceModel::stampHistory(const std::string &str) {

  std::string history = pism_username_prefix(m_grid->com) + (str + "\n");

  m_output_global_attributes.set_string("history",
                               history + m_output_global_attributes.get_string("history"));

}

void IceModel::check_minimum_ice_thickness() const {

  IceModelVec::AccessList list(m_ice_thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_ice_thickness(i, j) < 0.0) {
        throw RuntimeError::formatted("Thickness is negative at point i=%d, j=%d", i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

//! Check if the thickness of the ice is too large and extend the grid if necessary.
/*!
  Extends the grid such that the new one has 2 (two) levels above the ice.
 */
void IceModel::check_maximum_ice_thickness() const {
  const double Lz = m_grid->Lz();

  IceModelVec::AccessList list(m_ice_thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_ice_thickness(i, j) > Lz) {
        throw RuntimeError::formatted("Ice thickness (%7.4f m) exceeds the height"
                                      " of the computational box (%7.4f m) at i=%d, j=%d.",
                                      m_ice_thickness(i, j), Lz, i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
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
