// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

PetscErrorCode IceModel::computeMaxDiffusivity(bool updateDiffusViewer) {
  // NOTE:  Assumes IceModel::vubar, vvbar holds correct 
  // and up-to-date values of velocities

  PetscErrorCode ierr;

  PetscScalar **h, **H, **ubar, **vbar, **ub, **vb, **D, **mask;
  PetscScalar Dmax = 0.0;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        if ( (intMask(mask[i][j]) == MASK_SHEET) || 
             ((doSuperpose == PETSC_TRUE) && (modMask(mask[i][j]) != MASK_FLOATING)) ) {
          const PetscScalar h_x=(h[i+1][j]-h[i-1][j])/(2.0*grid.p->dx);
          const PetscScalar h_y=(h[i][j+1]-h[i][j-1])/(2.0*grid.p->dy);
          const PetscScalar alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
          // note: when basal sliding is proportional to surface slope, as
          // it usually will be when sliding occurs in a MASK_SHEET area, then
          //    D = H Ubar / alpha
          // is the correct formula; note division by zero is avoided by
          // addition to alpha
          PetscScalar Ubarmag;
          if (doSuperpose == PETSC_TRUE) {
            Ubarmag = sqrt(PetscSqr(ubar[i][j] - ub[i][j]) + PetscSqr(vbar[i][j] - vb[i][j]));
          } else {
            Ubarmag = sqrt(PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]));
          }
          const PetscScalar d =
               H[i][j] * Ubarmag/(alpha + DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE);
          if (d > Dmax) Dmax = d;
          D[i][j] = d;
        } else {
          D[i][j] = 0.0; // no diffusivity in non-SIA regions
        }
      } else {
        D[i][j] = 0.0; // no diffusivity if no ice
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);

  if (updateDiffusViewer) { // view diffusivity (m^2/s)
    ierr = update2DViewer('D',vWork2d[0],1.0); CHKERRQ(ierr);
  }

  ierr = PetscGlobalMax(&Dmax, &gDmax, grid.com); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::computeBasalDrivingStress(Vec myVec) {
  // puts   f_basal = \rho g H |\grad h|   into myVec; does not communicate ghosts
  PetscErrorCode ierr;

  PetscScalar **h, **H, **fbasal;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myVec, &fbasal); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar h_x = (h[i+1][j]-h[i-1][j])/(2.0*grid.p->dx);
      const PetscScalar h_y = (h[i][j+1]-h[i][j-1])/(2.0*grid.p->dy);
      const PetscScalar alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
      fbasal[i][j] = ice.rho * grav * H[i][j] * alpha;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, myVec, &fbasal); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::adaptTimeStepDiffusivity() {
  // note computeMaxDiffusivity() must be called before this to set gDmax
  // note that adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so 
  // dt <= adapt_ratio * dx^2/maxD (if dx=dy)
  // reference: Morton & Mayers 2nd ed. pp 62--63
  const PetscScalar  
          gridfactor = 1.0/(grid.p->dx*grid.p->dx) + 1.0/(grid.p->dy*grid.p->dy);
  dt_from_diffus = adaptTimeStepRatio
                     * 2 / ((gDmax + DEFAULT_ADDED_TO_GDMAX_ADAPT) * gridfactor);
  if ((doTempSkip == PETSC_TRUE) && (tempskipCountDown == 0)) {
    const PetscScalar  conservativeFactor = 0.8;
    // typically "dt" in next line is from CFL, but might be from other, e.g. maxdt
    tempskipCountDown = (PetscInt) floor(conservativeFactor * (dt / dt_from_diffus));
    tempskipCountDown = (tempskipCountDown > tempskipMax) ? tempskipMax : tempskipCountDown;
  } // if tempskipCountDown > 0 then it will get decremented at the mass balance step
  if (dt_from_diffus < dt) {
    dt = dt_from_diffus;
    adaptReasonFlag = 'd';
  }
  return 0;
}


PetscErrorCode IceModel::determineTimeStep(const bool doTemperatureCFL) {
  PetscErrorCode ierr;

  if ( (runtimeViewers[cIndex('D')] != PETSC_NULL) 
       || ( (doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE) ) ) {
    ierr = computeMaxDiffusivity(true); CHKERRQ(ierr);
  }
  const PetscScalar timeToEnd = (endYear-grid.p->year) * secpera;
  if (dt_force > 0.0) {
    dt = dt_force; // override usual dt mechanism
    adaptReasonFlag = 'f';
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
  } else {
    dt = maxdt;
    adaptReasonFlag = 'm';
    if ((doAdaptTimeStep == PETSC_TRUE) && (doTemp == PETSC_TRUE)
        && doTemperatureCFL) {
      // CFLmaxdt is set by computeMax3DVelocities() in call to velocity() iMvelocity.cc
      dt_from_cfl = CFLmaxdt;
      if (dt_from_cfl < dt) {
        dt = dt_from_cfl;
        adaptReasonFlag = 'c';
      }
    } 
    if ((doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE)
        && (useSSAVelocity)) {
      // CFLmaxdt2D is set by broadcastSSAVelocity()
      if (CFLmaxdt2D < dt) {
        dt = CFLmaxdt2D;
        adaptReasonFlag = 'u';
      }
    }
    if ((doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE)) {
      // note: if doTempSkip then tempskipCountDown = floor(dt_from_cfl/dt_from_diffus)
      ierr = adaptTimeStepDiffusivity(); CHKERRQ(ierr); // might set adaptReasonFlag = 'd'
    }
    if ((maxdt_temporary > 0.0) && (maxdt_temporary < dt)) {
      dt = maxdt_temporary;
      adaptReasonFlag = 't';
    }
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
    if ((adaptReasonFlag == 'm') || (adaptReasonFlag == 't') || (adaptReasonFlag == 'e')) {
      if (tempskipCountDown > 1) tempskipCountDown = 1; 
    }
  }    

/*
verbPrintf(1,grid.com,
   "\n   leaving determineTimeStep();\n     endYear, grid.p->year, dt_force = %f,%f,%e\n",
   endYear, grid.p->year, dt_force);
verbPrintf(1,grid.com,
   "     maxdt, CFLmaxdt, maxdt_temporary, gDmax = %e,%e,%e,%e\n",
   maxdt, CFLmaxdt, maxdt_temporary, gDmax);
verbPrintf(1,grid.com,
   "     dt, dt/secpera = %e, %f\n",
   dt, dt/secpera);
*/
  return 0;
}


PetscErrorCode IceModel::additionalAtStartTimestep() {
  // does nothing; allows derived classes to do more per-step computation,
  // reporting and checking, etc.
  return 0;
}


PetscErrorCode IceModel::additionalAtEndTimestep() {
  // does nothing; allows derived classes to do more per-step computation,
  // reporting and checking, etc.
  return 0;
}


PetscErrorCode IceModel::getHorSliceOf3D(Vec v3D, Vec &gslice, PetscInt k) {
  PetscErrorCode ierr;
  PetscScalar    ***val, **slice_val;

  ierr = DAVecGetArray(grid.da2, gslice, &slice_val); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      slice_val[i][j] = val[i][j][k];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, gslice, &slice_val); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::getSurfaceValuesOf3D(Vec v3D, Vec &gsurf) {
  PetscErrorCode ierr;
  PetscScalar    ***val, **H, **surf_val;

  ierr = DAVecGetArray(grid.da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  const PetscScalar dz = grid.p->dz;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/dz));
      if (ks == grid.p->Mz-1) {  // if at top of grid
        surf_val[i][j] = val[i][j][ks];
      } else {  // general case is to linearly interpolate to z=H[i][j]
        surf_val[i][j] = ((H[i][j] - ks * dz) * val[i][j][ks]
                         + ((ks+1) * dz - H[i][j]) * val[i][j][ks+1]) / dz;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::initFromOptions() {
  PetscErrorCode ierr;

  ierr = initFromOptions(PETSC_TRUE); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::initFromOptions(PetscTruth doHook) {
  PetscErrorCode ierr;
  PetscTruth inFileSet, bootstrapSet, bootstrapSetLegacy;
  char inFile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootstrapSet);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif_legacy", inFile,
                               PETSC_MAX_PATH_LEN, &bootstrapSetLegacy);
  CHKERRQ(ierr);
  
  if (bootstrapSet == PETSC_TRUE) {
    ierr = bootstrapFromFile_netCDF(inFile); CHKERRQ(ierr);
  } else if(bootstrapSetLegacy == PETSC_TRUE) {
    ierr = bootstrapFromFile_netCDF_legacyAnt(inFile); CHKERRQ(ierr);
  } else if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
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
    ierr = grid.rescale(my_Lx*1000.0, grid.p->Ly, grid.p->Lz); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Ly", &my_Ly, &LySet); CHKERRQ(ierr);
  if (LySet == PETSC_TRUE) {
    ierr = grid.rescale(grid.p->Lx, my_Ly*1000.0, grid.p->Lz); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lz", &my_Lz, &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_TRUE) {
    ierr = grid.rescale(grid.p->Lx, grid.p->Ly, my_Lz); CHKERRQ(ierr);
  }

  // make sure the first vertical velocities do not use junk from uninitialized basal melt rate.
  ierr = VecSet(vbasalMeltRate, 0.0); CHKERRQ(ierr);

  ierr = initBasalTillModel(); CHKERRQ(ierr);
    
  ierr = bedDefSetup(); CHKERRQ(ierr);

  ierr = initPDDFromOptions(); CHKERRQ(ierr);

  tempskipCountDown = 0;

  if (doHook == PETSC_TRUE) {
    ierr = afterInitHook(); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::afterInitHook() {
  PetscErrorCode ierr;
  PetscTruth     regridFileSet = PETSC_FALSE;
  char           regridFile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid", regridFile, PETSC_MAX_PATH_LEN,
                               &regridFileSet); CHKERRQ(ierr);
  if (regridFileSet == PETSC_TRUE) {
    ierr = regrid(regridFile); CHKERRQ(ierr);
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
  }

  ierr = stampHistoryCommand(); CHKERRQ(ierr);

  // initialization should be done; report on starting state
  ierr = verbPrintf(2,grid.com, 
           "  [computational box for ice: (%8.2f km) x (%8.2f km)",
           2*grid.p->Lx/1000.0,2*grid.p->Ly/1000.0); CHKERRQ(ierr);
  if (grid.p->Mbz > 1) {
    ierr = verbPrintf(2,grid.com,
         "\n                                 x (%8.2f m + %7.2f m bedrock)]\n"
         ,grid.p->Lz,grid.p->Lbz); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com," x (%8.2f m)]\n",grid.p->Lz); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,grid.com, 
           "  [grid cell dimensions     : (%8.2f km) x (%8.2f km) x (%8.2f m)]\n",
           grid.p->dx/1000.0,grid.p->dy/1000.0,grid.p->dz); CHKERRQ(ierr);

  // if -verbose then actually list all of IceParam
  ierr = verbPrintf(3,grid.com,"  IceParam: Mx = %d, My = %d, Mz = %d, Mbz = %d,\n",
                    grid.p->Mx,grid.p->My,grid.p->Mz,grid.p->Mbz); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "            Lx = %6.2f km, Ly = %6.2f m, Lz = %6.2f m, Lbz = %6.2f m,\n",
           grid.p->Lx/1000.0,grid.p->Ly/1000.0,grid.p->Lz,grid.p->Lbz); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "            dx = %6.3f km, dy = %6.3f km, dz = %6.3f m, year = %8.4f,\n",
           grid.p->dx/1000.0,grid.p->dy/1000.0,grid.p->dz,grid.p->year); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
     "            history = ****************\n%s            **************************\n"
     ,grid.p->history); CHKERRQ(ierr);
  
  ierr = createViewers(); CHKERRQ(ierr);
  return 0;
}


int IceModel::endOfTimeStepHook() {

  // SIGTERM makes PISM end, saving state under original -o name
  if (pism_signal == SIGTERM) {
    verbPrintf(1, grid.com, "Caught signal SIGTERM: exiting early.\n");
    return 1;
  }
  
  // SIGUSR1 makes PISM save state under filename based on the current year
  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "pism-%f.nc", grid.p->year);
    verbPrintf(1, grid.com, "Caught signal SIGUSR1: Writing intermediate file `%s'.\n",
               file_name);
    pism_signal = 0;
    dumpToFile_netCDF(file_name);
  }

  return 0;
}


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
  // ARGV.unshift($0).join(' ') and just ARGV.join(' ') if the executable name
  // was not needed.

  ierr = stampHistory(str); CHKERRQ(ierr);

  return 0;
}


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
  int length = snprintf(str, sizeof(str), "%s@%s %s : %s\n",
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


PetscErrorCode  IceModel::stampHistoryString(const char* string) {
  PetscErrorCode ierr;

  PetscInt historyLength = strlen(grid.p->history);
  PetscInt stringLength = strlen(string);

  if (stringLength + historyLength > (int)sizeof(grid.p->history) - 1) {
    ierr = PetscPrintf(grid.com, "Warning: History string overflow.  Truncating history.\n");
    CHKERRQ(ierr);

    // Don't overflow the buffer and null terminate.
    strncpy(grid.p->history, string, sizeof(grid.p->history));
    grid.p->history[sizeof(grid.p->history) - 1] = '\0';
  } else { // We are safe, so we can just write it.
    strcat(grid.p->history, string);
  }
  
  return 0;
}
