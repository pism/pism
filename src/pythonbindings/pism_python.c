// Copyright (C) 2011 David Maxwell
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

#include "pism_python.h"
#include "stdio.h"

PetscErrorCode globalMax(PetscReal local_max, PetscReal *result, MPI_Comm comm)
{
  return PetscGlobalMax(&local_max,result,comm);
}
PetscErrorCode globalMin(PetscReal local_min, PetscReal *result, MPI_Comm comm)
{
  return PetscGlobalMin(&local_min,result,comm);  
}
PetscErrorCode globalSum(PetscReal local_sum, PetscReal *result, MPI_Comm comm)
{
  return PetscGlobalSum(&local_sum,result,comm);  
}
