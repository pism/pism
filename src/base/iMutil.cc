// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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

#include <cstring>
#include <ctime>
#include "iceModel.hh"
#include "pism_signal.h"
#include <petscvec.h>


//! Compute scalar basal driving stress for diagnostic purposes.
/*!
Computes the scalar value of the driving stress at the base of the ice:
   \f[ f_{basal} = \rho g H |\nabla h| \f]
Saves it in \c myVec, which is treated as global.  (Does not 
communicate ghosts.)
 */
PetscErrorCode IceModel::computeBasalDrivingStress(Vec myVec) {
  PetscErrorCode ierr;

  PetscScalar **h, **H, **fbasal;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myVec, &fbasal); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar h_x = (h[i+1][j]-h[i-1][j])/(2.0*grid.dx);
      const PetscScalar h_y = (h[i][j+1]-h[i][j-1])/(2.0*grid.dy);
      const PetscScalar alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
      fbasal[i][j] = ice->rho * grav * H[i][j] * alpha;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, myVec, &fbasal); CHKERRQ(ierr);
  return 0;
}


//! Does nothing in \c IceModel, but allows derived classes to do more per-step computation.
PetscErrorCode IceModel::additionalAtStartTimestep() {
  return 0;
}


//! Does nothing in \c IceModel, but allows derived classes to do more per-step computation.
PetscErrorCode IceModel::additionalAtEndTimestep() {
  return 0;
}


//! Manages the initialization of IceModel, especially from input file options.
PetscErrorCode IceModel::initFromOptions() {
  PetscErrorCode ierr;

  ierr = initFromOptions(PETSC_TRUE); CHKERRQ(ierr);
  return 0;
}


//! Version of initFromOptions() which allows turning off afterInitHook().
PetscErrorCode IceModel::initFromOptions(PetscTruth doHook) {
  PetscErrorCode ierr;
  PetscTruth inFileSet, bootstrapSet; // bootstrapSetLegacy;
  char inFile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet);  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootstrapSet);  CHKERRQ(ierr);
  
  if ((inFileSet == PETSC_TRUE) && (bootstrapSet == PETSC_TRUE)) {
    SETERRQ(2,"both options '-if' and '-bif' set; not allowed");
  }

  if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile_netCDF(inFile); CHKERRQ(ierr);
  } else if (bootstrapSet == PETSC_TRUE) {
    ierr = bootstrapFromFile_netCDF(inFile); CHKERRQ(ierr);
  }

  if (! isInitialized()) {
    SETERRQ(1,"Model has not been initialized from a file or by a derived class.");
  }
  
  if (yearsStartRunEndDetermined == PETSC_FALSE) {
    ierr = setStartRunEndYearsFromOptions(PETSC_FALSE);  CHKERRQ(ierr);
  }
  
  // runtime options take precedence in setting of -Lx,-Ly,-Lz *including*
  // if initialization is from an input file
  PetscTruth     LxSet, LySet, LzSet;
  PetscScalar    my_Lx, my_Ly, my_Lz;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lx", &my_Lx, &LxSet); CHKERRQ(ierr);
  if (LxSet == PETSC_TRUE) {
    ierr = grid.rescale_using_zlevels(my_Lx*1000.0, grid.Ly); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Ly", &my_Ly, &LySet); CHKERRQ(ierr);
  if (LySet == PETSC_TRUE) {
    ierr = grid.rescale_using_zlevels(grid.Lx, my_Ly*1000.0); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lz", &my_Lz, &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "resetting vertical levels base on options and user option -Lz ...\n");
      CHKERRQ(ierr);
    ierr = determineSpacingTypeFromOptions(); CHKERRQ(ierr);
    ierr = grid.rescale_and_set_zlevels(grid.Lx, grid.Ly, my_Lz); CHKERRQ(ierr);
  }

  // make sure the first vertical velocities do not use junk from 
  //   uninitialized basal melt rate.
  ierr = VecSet(vbasalMeltRate, 0.0); CHKERRQ(ierr);

  ierr = initBasalTillModel(); CHKERRQ(ierr);
    
  ierr = bedDefSetup(); CHKERRQ(ierr);

  ierr = initPDDFromOptions(); CHKERRQ(ierr);

  ierr = initForcingFromOptions(); CHKERRQ(ierr);

  tempskipCountDown = 0;

  if (doHook == PETSC_TRUE) {
    ierr = afterInitHook(); CHKERRQ(ierr);
  }

  return 0;
}


//! Complete initialization: regrid if desired, report on computational domain and grid, create viewers.
PetscErrorCode IceModel::afterInitHook() {
  PetscErrorCode ierr;

  PetscTruth     regridFileSet = PETSC_FALSE;
  char           regridFile[PETSC_MAX_PATH_LEN];

  // initialization should be done by here!
  
  // report on computational box
  ierr = verbPrintf(2,grid.com, 
           "  [computational box for ice: (%8.2f km) x (%8.2f km)",
           2*grid.Lx/1000.0,2*grid.Ly/1000.0); CHKERRQ(ierr);
  if (grid.Mbz > 1) {
    ierr = verbPrintf(2,grid.com,
         "\n                                 x (%8.2f m + %7.2f m bedrock)]\n"
         ,grid.Lz,grid.Lbz); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com," x (%8.2f m)]\n",grid.Lz); CHKERRQ(ierr);
  }
  
  // report on grid cell dims
  if (grid.isEqualVertSpacing()) {
    ierr = verbPrintf(2,grid.com, 
           "  [grid cell dims (equal dz): (%8.2f km) x (%8.2f km) x (%8.2f m)]\n",
           grid.dx/1000.0,grid.dy/1000.0,grid.dzMIN); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, 
           "  [hor. grid cell dimensions: (%8.2f km) x (%8.2f km)]\n",
           grid.dx/1000.0,grid.dy/1000.0); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, 
           "  [vertical grid spacing in ice not equal: %.3f m < dz < %.3f m]\n",
           grid.dzMIN,grid.dzMAX); CHKERRQ(ierr);
    PetscInt    myMz, dummyM;
    ierr = getMzMbzForTempAge(myMz,dummyM); CHKERRQ(ierr);
    if (myMz > 1000) {
      ierr = verbPrintf(1,grid.com,
        "\n\n WARNING: Using more than 1000 vertical levels internally\n"
        "   in temperatureStep()!\n\n");  CHKERRQ(ierr);
    }
    ierr = verbPrintf(2,grid.com, 
           "  [fine equal spacing used in temperatureStep(): Mz = %d, dzEQ = %.3f m]\n",
           myMz,grid.Lz / ((PetscScalar) (myMz - 1))); CHKERRQ(ierr);
    if (grid.Mbz > 1) {
      ierr = verbPrintf(2,grid.com, 
         "  [vertical spacing in bedrock: dz = %.3f m for bottom layer\n"
         "     and dz = %.3f m for top layer]\n",
         grid.zblevels[1]-grid.zblevels[0],
         grid.zblevels[grid.Mbz-1]-grid.zblevels[grid.Mbz-2]); CHKERRQ(ierr);
    }
  }

  // if -verbose then actually list all members of grid
  ierr = verbPrintf(3,grid.com,
           "            Mx = %d, My = %d, Mz = %d, Mbz = %d,\n",
                    grid.Mx,grid.My,grid.Mz,grid.Mbz); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "            Lx = %6.2f km, Ly = %6.2f m, Lz = %6.2f m, Lbz = %6.2f m,\n",
           grid.Lx/1000.0,grid.Ly/1000.0,grid.Lz,grid.Lbz); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "            dx = %6.3f km, dy = %6.3f km, year = %8.4f,\n",
           grid.dx/1000.0,grid.dy/1000.0,grid.year); CHKERRQ(ierr);
  ierr = grid.printVertLevels(5); CHKERRQ(ierr);  // only if verbose 5
  ierr = verbPrintf(5,grid.com,
     "            history = ****************\n%s            **************************\n"
     ,history); CHKERRQ(ierr);
  
  // miscellaneous
  ierr = stampHistoryCommand(); CHKERRQ(ierr);
  ierr = createViewers(); CHKERRQ(ierr);

  // read new values from regrid file, and overwrite current, if desired
  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid", regridFile, PETSC_MAX_PATH_LEN,
                               &regridFileSet); CHKERRQ(ierr);
  if (regridFileSet == PETSC_TRUE) {
    ierr = regrid_netCDF(regridFile); CHKERRQ(ierr);
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
  }

  return 0;
}


//! Catch signals -USR1 and -TERM; in the former case save and continue; in the latter, save and stop.
/*!
Signal \c SIGTERM makes PISM end, saving state under original \c -o name 
(or default name).  We also add an indication to the history attribute 
of the output NetCDF file.

Signal \c SIGUSR1 makes PISM save state under a filename based on the
the name of the executable (e.g. \c pismr or \c pismv) and the current 
model year.  There is no indication in the history attribute of the output 
NetCDF file because there is no effect.  There is an indication at \c stdout.
 */
int IceModel::endOfTimeStepHook() {
  
  if (pism_signal == SIGTERM) {
    verbPrintf(1, grid.com, 
       "Caught signal SIGTERM:  EXITING EARLY and saving with original filename.\n");
    char str[HISTORY_STRING_LENGTH];
    snprintf(str, sizeof(str), 
       "EARLY EXIT caused by signal SIGTERM.  Completed timestep at year=%.3f.",
       grid.year);
    stampHistory(str);
    return 1;
  }
  
  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "%s-%5.3f.nc",
             executable_short_name, grid.year);
    verbPrintf(1, grid.com, 
       "Caught signal SIGUSR1:  Writing intermediate file `%s'.\n",
       file_name);
    pism_signal = 0;
    dumpToFile_netCDF(file_name);
  }

  return 0;
}


//! Build a history string from the command which invoked PISM.
PetscErrorCode  IceModel::stampHistoryCommand() {
  PetscErrorCode ierr;
  PetscInt argc;
  char **argv;
  
  ierr = PetscGetArgs(&argc, &argv); CHKERRQ(ierr);
  
  char str[HISTORY_STRING_LENGTH];
  
  strncpy(str, argv[0], sizeof(str)); // Does not null terminate on overflow
  str[sizeof(str) - 1] = '\0';
  for (PetscInt i=1; i < argc; i++) {
    PetscInt remaining_bytes = sizeof(str) - strlen(str) - 1;
    // strncat promises to null terminate, so we must only make sure that the
    // end of the buffer is not overwritten.
    strncat(str, " ", remaining_bytes--);
    strncat(str, argv[i], remaining_bytes);
  }

  // All this hooplah is the C equivalent of the Ruby expression
  // ARGV.unshift($0).join(' ') and just ARGV.join(' ')
  // if the executable name was not needed.

  ierr = stampHistory(str); CHKERRQ(ierr);

  return 0;
}


//! Build the particular history string associated to the end of a PISM run.
PetscErrorCode  IceModel::stampHistoryEnd() {
  PetscErrorCode ierr;
  PetscLogDouble flops, my_flops;
  char str[HISTORY_STRING_LENGTH];
  MPI_Datatype mpi_type;

  ierr = PetscGetFlops(&my_flops); CHKERRQ(ierr);

  ierr = PetscDataTypeToMPIDataType(PETSC_DOUBLE, &mpi_type); CHKERRQ(ierr);
  MPI_Reduce(&my_flops, &flops, 1, mpi_type, MPI_SUM, 0, grid.com);
  
  snprintf(str, sizeof(str), "Done. Petsc reports %.2f MFlops on %d processes.",
           flops * 1.0e-6, grid.size);

  ierr = stampHistory(str); CHKERRQ(ierr);
  
  return 0;
}


//! Get time and user/host name and add it to the given string.  Then call stampHistoryString(). 
PetscErrorCode  IceModel::stampHistory(const char* string) {
  PetscErrorCode ierr;

  time_t now;
  tm tm_now;
  now = time(NULL);
  localtime_r(&now, &tm_now);

  char date_str[50];
  // Format specifiers for strftime():
  //   %F : ISO date format
  //   %T : Full 24 hour time
  //   %Z : Time Zone name
  //   %z : Time zone offset
  strftime(date_str, sizeof(date_str), "%F %T %Z", &tm_now);

  char username[50];
  ierr = PetscGetUserName(username, sizeof(username)); CHKERRQ(ierr);
  char hostname[100];
  ierr = PetscGetHostName(hostname, sizeof(hostname)); CHKERRQ(ierr);
  
  char str[HISTORY_STRING_LENGTH];
  int length = snprintf(str, sizeof(str), "%s@%s %s:  %s\n",
                        username, hostname, date_str, string);
  
  if (length < 0) {
    SETERRQ(1, "Output error or snprintf() is not C99 compliant.");
    // Old implementations of snprintf() will return `-1' when the string is
    // truncated.  If this is the case on some platform, we need to change this
    // check to allow for that possibility.
  }
  if (length > (int)sizeof(str)) {
    ierr = PetscPrintf(grid.com,
       "Warning: command line truncated by %d chars in history.\n",
       length + 1 - sizeof(str)); CHKERRQ(ierr);
    str[sizeof(str) - 2] = '\n';
    str[sizeof(str) - 1] = '\0';
  }

  ierr = stampHistoryString(str); CHKERRQ(ierr);
  
  return 0;
}


//! Add the given string to the history data member in IceModel.
PetscErrorCode  IceModel::stampHistoryString(const char* string) {
  PetscErrorCode ierr;

  PetscInt historyLength = strlen(history);
  PetscInt stringLength = strlen(string);

  if (stringLength + historyLength > (int)sizeof(history) - 1) {
    ierr = PetscPrintf(grid.com, 
         "Warning: History string overflow.  Truncating history.\n");
         CHKERRQ(ierr);

    // Don't overflow the buffer and null terminate.
    strncpy(history, string, sizeof(history));
    history[sizeof(history) - 1] = '\0';
  } else { // We are safe, so we can just write it.
    //OLD METHOD: append the latest command:    
    // strcat(history, string);
    //NEW METHOD: prepend it; this matches NCO behavior so commands are in order
    char tempstr[HISTORY_STRING_LENGTH];
    strcpy(tempstr,string);
    strcat(tempstr,history);
    strcpy(history,tempstr);
  }
  
  return 0;
}
