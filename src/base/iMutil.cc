// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <petscvec.h>

#include "iceModel.hh"
#include "pism_signal.h"
#include "PISMSurface.hh"
#include "PISMStressBalance.hh"
#include "enthalpyConverter.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "PISMDiagnostic.hh"
#include "bedrockThermalUnit.hh"

namespace pism {


//! Virtual.  Does nothing in `IceModel`.  Derived classes can do more computation in each time step.
PetscErrorCode IceModel::additionalAtStartTimestep() {
  return 0;
}


//! Virtual.  Does nothing in `IceModel`.  Derived classes can do more computation in each time step.
PetscErrorCode IceModel::additionalAtEndTimestep() {
  return 0;
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
  PetscErrorCode ierr;
  
  if (pism_signal == SIGTERM) {
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGTERM:  EXITING EARLY and saving with original filename.\n");
    char str[TEMPORARY_STRING_LENGTH];
    snprintf(str, sizeof(str), 
       "EARLY EXIT caused by signal SIGTERM.  Completed timestep at time=%s.",
             grid.time->date().c_str());
    stampHistory(str);
    return 1;
  }
  
  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "%s-%s.nc",
             executable_short_name.c_str(), grid.time->date().c_str());
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGUSR1:  Writing intermediate file `%s' and flushing time series.\n\n",
       file_name);
    pism_signal = 0;
    dumpToFile(file_name);

    // flush all the time-series buffers:
    ierr = flush_timeseries(); CHKERRQ(ierr); 
  }

  if (pism_signal == SIGUSR2) {
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGUSR2:  Flushing time series.\n\n");
    pism_signal = 0;

    // flush all the time-series buffers:
    ierr = flush_timeseries(); CHKERRQ(ierr);
  }

  return 0;
}


//! Build a history string from the command which invoked PISM.
PetscErrorCode  IceModel::stampHistoryCommand() {
  PetscErrorCode ierr;
  
  char startstr[TEMPORARY_STRING_LENGTH];

  snprintf(startstr, sizeof(startstr), 
           "PISM (%s) started on %d procs.", PISM_Revision, (int)grid.size);
  ierr = stampHistory(std::string(startstr)); CHKERRQ(ierr);

  global_attributes.set_string("history",
                               pism_args_string() + global_attributes.get_string("history"));

  return 0;
}

PetscErrorCode IceModel::update_run_stats() {
  PetscErrorCode ierr;

  MPI_Datatype mpi_type;
  ierr = PetscDataTypeToMPIDataType(PETSC_DOUBLE, &mpi_type); CHKERRQ(ierr);

  // timing stats
  PetscLogDouble current_time, my_current_time;
  double wall_clock_hours, proc_hours, mypph;
  ierr = GetTime(&my_current_time); CHKERRQ(ierr);
  MPI_Allreduce(&my_current_time, &current_time, 1, mpi_type, MPI_MAX, grid.com);

  wall_clock_hours = (current_time - start_time) / 3600.0;

  proc_hours = grid.size * wall_clock_hours;

  // MYPPH stands for "model years per processor hour"
  mypph = grid.convert(grid.time->current() - grid.time->start(), "seconds", "years") / proc_hours;

  MPI_Bcast(&mypph, 1, MPI_DOUBLE, 0, grid.com);

  // get PETSc's reported number of floating point ops (*not* per time) on this
  //   process, then sum over all processes
  PetscLogDouble flops, my_flops;
  ierr = PetscGetFlops(&my_flops); CHKERRQ(ierr);
  MPI_Allreduce(&my_flops, &flops, 1, mpi_type, MPI_SUM, grid.com);

  run_stats.set_double("wall_clock_hours", wall_clock_hours);
  run_stats.set_double("processor_hours", proc_hours);
  run_stats.set_double("model_years_per_processor_hour", mypph);
  run_stats.set_double("PETSc_MFlops", flops * 1.0e-6);
  run_stats.set_double("grid_dx_meters", grid.dx);
  run_stats.set_double("grid_dy_meters", grid.dy);
  run_stats.set_double("grid_dz_min_meters", grid.dzMIN);
  run_stats.set_double("grid_dz_max_meters", grid.dzMAX);
  if (btu != NULL) {
    run_stats.set_double("grid_dzb_meters", btu->get_vertical_spacing());
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

  return 0;
}

//! Build the particular history string associated to the end of a PISM run,
//! including a minimal performance assessment.
PetscErrorCode  IceModel::stampHistoryEnd() {
  PetscErrorCode ierr;

  ierr = update_run_stats(); CHKERRQ(ierr);

  // build and put string into global attribute "history"
  char str[TEMPORARY_STRING_LENGTH];

  snprintf(str, TEMPORARY_STRING_LENGTH,
    "PISM done.  Performance stats: %.4f wall clock hours, %.4f proc.-hours, %.4f model years per proc.-hour, PETSc MFlops = %.2f.",
           run_stats.get_double("wall_clock_hours"),
           run_stats.get_double("processor_hours"),
           run_stats.get_double("model_years_per_processor_hour"),
           run_stats.get_double("PETSc_MFlops"));

  ierr = stampHistory(str); CHKERRQ(ierr);

  return 0;
}


//! Get time and user/host name and add it to the given string.
PetscErrorCode  IceModel::stampHistory(const std::string &str) {

  std::string history = pism_username_prefix(grid.com) + (str + "\n");

  global_attributes.set_string("history",
                               history + global_attributes.get_string("history"));
  
  return 0;
}

//! Check if the thickness of the ice is too large and extend the grid if necessary.
/*!
  Extends the grid such that the new one has 2 (two) levels above the ice.
 */
PetscErrorCode IceModel::check_maximum_thickness() {
  PetscErrorCode  ierr;
  double H_min, H_max;

  ierr = ice_thickness.range(H_min, H_max); CHKERRQ(ierr);
  if (grid.Lz >= H_max) return 0;

  ierr = PetscPrintf(grid.com,
                     "PISM ERROR: Max ice thickness (%7.4f m) is greater than the height of the computational box (%7.4f m).\n",
                     H_max, grid.Lz); CHKERRQ(ierr);
  PISMEnd();

  return 0;
}


//! Allows derived classes to extend their own IceModelVec3's in vertical.
/*! Base class version does absolutely nothing. */
PetscErrorCode IceModel::check_maximum_thickness_hook(const int /*old_Mz*/) {
  return 0;
}

} // end of namespace pism
