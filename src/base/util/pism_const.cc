// Copyright (C) 2007--2017 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include <petscfix.h>
#include <petsctime.h>
#include <petscsys.h>

#include <sstream>
#include <ctime>
#include <algorithm>
#include <sys/types.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>

#include "pism_const.hh"
#include "pism_utilities.hh"
#include "error_handling.hh"

namespace pism {

PetscLogDouble GetTime() {
  PetscLogDouble result;
  PetscErrorCode ierr = PetscTime(&result); PISM_CHK(ierr, "PetscTime");
  return result;
}

} // end of namespace pism
