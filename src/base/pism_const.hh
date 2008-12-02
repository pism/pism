// Copyright (C) 2007--2008 Jed Brown, Ed Bueler and Constantine Khroulev
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

const PetscScalar gasConst_R = 8.31441;      // J/(mol K)    Gas Constant
const PetscScalar grav       = 9.81;         // m/s^2        acceleration of gravity
const PetscScalar secpera    = 3.1556926e7;
const PetscScalar pi         = 3.14159265358979;

const PetscInt TEMPORARY_STRING_LENGTH = 32768; // 32KiB ought to be enough.

PetscErrorCode getFlowLawNumber(PetscInt &flowLawNum, const PetscInt defaultFLN);
PetscErrorCode userChoosesIceType(MPI_Comm com, IceType* &ice);
PetscErrorCode userChoosesIceType(MPI_Comm com, IceType* &ice, const PetscInt defaultFLN);

PetscErrorCode setVerbosityLevel(PetscInt level);
PetscInt getVerbosityLevel();
PetscErrorCode verbosityLevelFromOptions();
PetscErrorCode verbPrintf(const int thresh, MPI_Comm comm,const char format[],...);
int compare_doubles (const void *a, const void *b);	// for sorting

#endif
