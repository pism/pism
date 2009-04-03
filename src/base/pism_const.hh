// Copyright (C) 2007--2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __pism_const_hh
#define __pism_const_hh

#include <petsc.h>
#include "materials.hh"

extern const char *PISM_Revision;

const PetscScalar gasConst_R = 8.31441;      // J/(mol K)    Gas Constant
const PetscScalar earth_grav = 9.81;         // m/s^2        acceleration of gravity
const PetscScalar secpera    = 3.15569259747e7; // The constant used in UDUNITS
						// (src/udunits/pismudunits.dat)
const PetscScalar pi         = 3.14159265358979;


// following numerical values have some significance; see updateSurfaceElevationAndMask()
enum PismMask {
  MASK_SHEET = 1,
  MASK_DRAGGING = 2,
  MASK_FLOATING = 3,
  // (PismModMask(mask[i][j]) == MASK_FLOATING) is criteria for floating; ..._OCEAN0 only used if -ocean_kill
  MASK_FLOATING_OCEAN0 = 7
};

// note no range checking in these two:
static inline int PismIntMask(PetscScalar maskvalue) {
  return static_cast<int>(floor(maskvalue + 0.5));
}

static inline int PismModMask(PetscScalar maskvalue) {
  int intmask = static_cast<int>(floor(maskvalue + 0.5));
  if (intmask > MASK_FLOATING) {
    return intmask - 4;
  } else {
    return intmask;
  }
}


// Standard C++ does not have a "NaN", or an "isnan()".  We need an alternative
//   and this is an admittedly lame one.  If there is a realiable way to check if
//   IEEE 754 NAN is available, then PISM_NAN could be set to that.  Note most
//   negative number in IEEE double is ~= -1.8e308, and that if PetscScalar !=
//   double, use of this may through an error (and that would probably be good).
const double PISM_DOUBLE_NAN = -1.234567890123456e308;

const PetscInt TEMPORARY_STRING_LENGTH = 32768; // 32KiB ought to be enough.

PetscErrorCode setVerbosityLevel(PetscInt level);
PetscInt getVerbosityLevel();
PetscErrorCode verbosityLevelFromOptions();

// handy functions for processing options:
PetscErrorCode check_option(const char name[], PetscTruth &flag);
PetscErrorCode ignore_option(MPI_Comm com, const char name[]);
PetscErrorCode check_old_option_and_stop(MPI_Comm com, const char old_name[], const char new_name[]);
PetscErrorCode stop_if_set(MPI_Comm com, const char name[]);

PetscErrorCode verbPrintf(const int thresh, MPI_Comm comm,const char format[],...);
int compare_doubles (const void *a, const void *b);	// for sorting

#endif
