// Copyright (C) 2012 PISM Authors
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

#ifndef _PISM_PETSC32_COMPAT_H_
#define _PISM_PETSC32_COMPAT_H_

#include <petsc.h>

// PETSc 3.2 compatibility header
//
// This takes care of some changes in the PETSc API between versions 3.2 and
// 3.3; it should be included in .cc files calling functions that were renamed.

#if !defined(PETSC_VERSION_LT)
#  define PISM_PETSC32_COMPAT 1
#elif PETSC_VERSION_LT(3,3,0)
#  define PISM_PETSC32_COMPAT 1
#else
#  define PISM_PETSC32_COMPAT 0
#endif

#if PISM_PETSC32_COMPAT==1
# define PetscObjectTypeCompare(obj,type,flag) PetscTypeCompare(obj,type,flag)
# define DMCreateMatrix(a,b,c) DMGetMatrix(a,b,c)
# define SNESDMComputeFunction SNESDAFormFunction
# define SNESDMComputeJacobian SNESDAComputeJacobian
#endif

#endif /* _PISM_PETSC32_COMPAT_H_ */

