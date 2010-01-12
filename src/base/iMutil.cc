// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "iceModel.hh"
#include "pism_signal.h"
#include <petscvec.h>

//! Virtual.  Does nothing in \c IceModel.  Derived classes can do more computation in each time step.
PetscErrorCode IceModel::additionalAtStartTimestep() {
  return 0;
}


//! Virtual.  Does nothing in \c IceModel.  Derived classes can do more computation in each time step.
PetscErrorCode IceModel::additionalAtEndTimestep() {
  return 0;
}

//! Catch signals -USR1 and -TERM; in the former case save and continue; in the latter, save and stop.
/*!
Signal \c SIGTERM makes PISM end, saving state under original \c -o name 
(or default name).  We also add an indication to the history attribute 
of the output NetCDF file.

Signal \c SIGUSR1 makes PISM save state under a filename based on the
the name of the executable (e.g. \c pismr or \c pismv) and the current 
model year.  In addition the time series (\c -ts_file, etc.) is flushed out
There is no indication of these actions in the history attribute of the output (\c -o)
NetCDF file because there is no effect on it, but there is an indication at \c stdout.
 */
int IceModel::endOfTimeStepHook() {
  
  if (pism_signal == SIGTERM) {
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGTERM:  EXITING EARLY and saving with original filename.\n");
    char str[TEMPORARY_STRING_LENGTH];
    snprintf(str, sizeof(str), 
       "EARLY EXIT caused by signal SIGTERM.  Completed timestep at year=%.3f.",
       grid.year);
    stampHistory(str);
    return 1;
  }
  
  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "%s-%5.3f.nc",
             executable_short_name.c_str(), grid.year);
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGUSR1:  Writing intermediate file `%s' and flushing time series.\n\n",
       file_name);
    pism_signal = 0;
    dumpToFile(file_name);

    // flush all the time-series buffers:
    vector<DiagnosticTimeseries*>::iterator i;
    for (i = timeseries.begin(); i < timeseries.end(); ++i) {
      (*i)->flush();
    }
  }

  return 0;
}


//! Build a history string from the command which invoked PISM.
PetscErrorCode  IceModel::stampHistoryCommand() {
  PetscErrorCode ierr;
  PetscInt argc;
  char **argv;
  
  ierr = PetscGetArgs(&argc, &argv); CHKERRQ(ierr);
  
  char startstr[TEMPORARY_STRING_LENGTH];

  snprintf(startstr, sizeof(startstr), 
           "PISM (%s) started on %d procs.", PISM_Revision, (int)grid.size);
  ierr = stampHistory(string(startstr)); CHKERRQ(ierr);

  // Create a string with space-separated command-line arguments:
  string cmdstr;
  for (int j = 0; j < argc; j++)
    cmdstr += string(" ") + argv[j];
  cmdstr += "\n";

  global_attributes.prepend_history(cmdstr);

  return 0;
}


//! Build the particular history string associated to the end of a PISM run.
PetscErrorCode  IceModel::stampHistoryEnd() {
  PetscErrorCode ierr;
  PetscLogDouble flops, my_flops;
  char str[TEMPORARY_STRING_LENGTH];
  MPI_Datatype mpi_type;

  ierr = PetscGetFlops(&my_flops); CHKERRQ(ierr);

  ierr = PetscDataTypeToMPIDataType(PETSC_DOUBLE, &mpi_type); CHKERRQ(ierr);
  MPI_Reduce(&my_flops, &flops, 1, mpi_type, MPI_SUM, 0, grid.com);
  
  snprintf(str, sizeof(str), "PISM done.  PETSc MFlops = %.2f.",
           flops * 1.0e-6);

  ierr = stampHistory(str); CHKERRQ(ierr);
  
  return 0;
}


//! Get time and user/host name and add it to the given string.
PetscErrorCode  IceModel::stampHistory(string str) {

  global_attributes.prepend_history(username_prefix() + (str + "\n"));
  
  return 0;
}

//! Check if the thickness of the ice is too large and extend the grid if necessary.
/*!
  Extends the grid such that the new one has 2 (two) levels above the ice.
 */
PetscErrorCode IceModel::check_maximum_thickness() {
  PetscErrorCode  ierr;
  PetscReal H_min, H_max, dz_top;
  double *new_zlevels, *new_zblevels;
  const int old_Mz = grid.Mz;
  int N = 0; 			// the number of new levels

  ierr = vH.range(H_min, H_max); CHKERRQ(ierr);
  if (grid.Lz >= H_max) return 0;

  if (grid.initial_Mz == 0)
    grid.initial_Mz = grid.Mz;
  else if (grid.Mz > grid.initial_Mz * 2) {
    ierr = PetscPrintf(grid.com,
		       "\n"
		       "PISM ERROR: Max ice thickness (%7.4f m) is greater than the height of the computational box (%7.4f m)"
		       " AND the grid has twice the initial number of vertical levels (%d) already. Exiting...\n",
		       H_max, grid.Lz, grid.initial_Mz); CHKERRQ(ierr);
    PetscEnd();
  }

  // So, we need to extend the grid. We find dz at the top of the grid,
  // create new zlevels and zblevels, then extend all the IceModelVec3s.

  // Find the vertical grid spacing at the top of the grid:
  dz_top = grid.Lz - grid.zlevels[old_Mz - 2];

  // Find the number of new levels:
  while (grid.Lz + N * dz_top <= H_max) N++;
  // This makes sure that we always have *exactly* two levels strictly above the ice:
  if (grid.Lz + N * dz_top > H_max) N += 1;
  else N += 2;

  ierr = verbPrintf(2, grid.com,
		    "\n"
		    "PISM WARNING: max ice thickness (%7.4f m) is greater than the height of the computational box (%7.4f m)...\n"
		    "              Adding %d new grid levels %7.4f m apart...\n",
		    H_max, grid.Lz, N, dz_top); CHKERRQ(ierr);

  // Create new zlevels and zblevels:
  new_zblevels = new double[grid.Mbz];
  new_zlevels = new double[old_Mz + N];

  for (int j = 0; j < grid.Mbz; j++)
    new_zblevels[j] = grid.zblevels[j];
  for (int j = 0; j < old_Mz; j++)
    new_zlevels[j] = grid.zlevels[j];

  // Fill the new levels:
  for (int j = 0; j < N; j++)
    new_zlevels[old_Mz + j] = grid.Lz + dz_top * (j + 1);

  ierr = grid.set_vertical_levels(old_Mz + N, grid.Mbz,
				  new_zlevels, new_zblevels); CHKERRQ(ierr);
  delete[] new_zlevels;
  delete[] new_zblevels;

  // Done with the grid. Now we need to extend IceModelVec3s.

  // We use surface temperatures to extend T3 and Tnew3. We get them from the
  // PISMAtmosphereCoupler.

  IceModelVec2 *pccTs;

  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, 0.0, pccTs); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }

  // Model state 3D vectors:
  ierr =     u3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     v3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     w3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr = Sigma3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     T3.extend_vertically(old_Mz, *pccTs); CHKERRQ(ierr);

  // Work 3D vectors:
  ierr =         Tnew3.extend_vertically(old_Mz, *pccTs); CHKERRQ(ierr);
  ierr = Sigmastag3[0].extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr = Sigmastag3[1].extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     Istag3[0].extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     Istag3[1].extend_vertically(old_Mz, 0); CHKERRQ(ierr);

  // deal with 3D age conditionally
  if (config.get_flag("do_age")) {
    ierr = tau3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
    ierr = taunew3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::report_grid_parameters() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com, "computational domain and grid:\n"); CHKERRQ(ierr);
  // report on computational box
  ierr = verbPrintf(2,grid.com, 
           "                    spatial domain   %.2f km x %.2f km",
           2*grid.Lx/1000.0,2*grid.Ly/1000.0); CHKERRQ(ierr);
  if (grid.Mbz > 1) {
    ierr = verbPrintf(2,grid.com," x (%.2f m + %.2f m bedrock)\n"
         ,grid.Lz,grid.Lbz); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com," x %.2f m\n",grid.Lz); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com,
           "                     time interval   start = %.2f a, end = %.2f a; run length = %.2f a\n",
		    grid.start_year, grid.end_year, grid.end_year - grid.start_year);
  
  // report on grid cell dims
  if (grid.ice_vertical_spacing == EQUAL) {
    ierr = verbPrintf(2,grid.com, 
           "  grid cell dims (equal dz in ice)   %.2f km x %.2f km x %.2f m\n",
           grid.dx/1000.0,grid.dy/1000.0,grid.dzMIN); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, 
           "   horizontal grid cell dimensions   %.2f km x %.2f km\n",
           grid.dx/1000.0,grid.dy/1000.0); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, 
           "      vertical grid spacing in ice   uneven, %d levels, %.3f m < dz < %.3f m\n",
		      grid.Mz, grid.dzMIN, grid.dzMAX); CHKERRQ(ierr);
  }

  PetscInt    fMz = 0;	// will be initialized by the call below
  if (grid.Mbz > 1) {
    if (grid.bed_vertical_spacing == EQUAL) {
      ierr = verbPrintf(2,grid.com, 
           "  vertical grid spacing in bedrock   equal, dz = %.3f m\n",
			grid.zblevels[1]-grid.zblevels[0]); CHKERRQ(ierr);
    } else {
    ierr = verbPrintf(2,grid.com, 
           "  vertical grid spacing in bedrock   uneven, %d levels, %.3f m < dz < %.3f m\n",
		      grid.Mbz, grid.dzbMIN, grid.dzbMAX); CHKERRQ(ierr);
    }
    PetscInt fMbz;
    PetscScalar fdz, fdzb, *fzlev, *fzblev;
    ierr = grid.get_fine_vertical_grid(fMz, fMbz, fdz, fdzb, fzlev, fzblev); CHKERRQ(ierr);
    delete[] fzlev; delete[] fzblev;
    ierr = verbPrintf(3,grid.com, 
           "   fine spacing used in energy/age   fMz = %d, fdz = %.3f m, fMbz = %d, fdzb = %.3f m\n",
           fMz, fdz, fMbz, fdzb); CHKERRQ(ierr);
  } else { // no bedrock case
    PetscScalar fdz, *fzlev;
    ierr = grid.get_fine_vertical_grid_ice(fMz, fdz, fzlev); CHKERRQ(ierr);
    delete[] fzlev;
    ierr = verbPrintf(3,grid.com, 
           "   fine spacing used in energy/age   fMz = %d, fdz = %.3f m\n",
           fMz, fdz); CHKERRQ(ierr);
  }
  if (fMz > 1000) {
    ierr = verbPrintf(1,grid.com,
      "\n\nWARNING: Using more than 1000 ice vertical levels internally in energy/age computation!\n\n");
      CHKERRQ(ierr);
  }

  // if -verbose (=-verbose 3) then (somewhat redundantly) list parameters of grid
  ierr = grid.printInfo(3); CHKERRQ(ierr);

  // if -verbose 5 then more stuff
  ierr = verbPrintf(5,grid.com,
       "  REALLY verbose output on IceGrid:\n"); CHKERRQ(ierr);
  ierr = grid.printVertLevels(5); CHKERRQ(ierr);
  
  return 0;
}

