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

#include "pism_python.hh"
#include "petsc.h"
#include "pism_python_signal.hh"
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

PetscErrorCode optionsGroupBegin(MPI_Comm comm,const char *prefix,const char *mess,const char *sec)
{
  PetscOptionsPublishCount=(PetscOptionsPublish?-1:1);
  return PetscOptionsBegin_Private(comm,prefix,mess,sec);
}

void optionsGroupNext()
{
  PetscOptionsPublishCount++;
}

bool optionsGroupContinue()
{
  return PetscOptionsPublishCount<2;
}

PetscErrorCode optionsGroupEnd()
{
  return PetscOptionsEnd_Private();
}


void set_abort_on_sigint(bool abort)
{
  gSIGINT_is_fatal = abort;
}

