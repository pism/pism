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

#ifndef __iceModelpreamble_hh
#define __iceModelpreamble_hh

// this header file is simply added to unclutter iceModel.hh

#include <petsc.h>
#include "materials.hh"

// remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond
 

struct titleNname {
  char title[100]; // these short titles appear on PETSc graphical viewers and in Matlab output file
  char name[30];   // these names are for Matlab output vars
};


struct PolarStereoParams {
  // these are "double" and not "float" ultimately because of how ncgen works
  double svlfp, // straight_vertical_longitude_from_pole; defaults to 0
         lopo,  // latitude_of_projection_origin; defaults to 90
         sp;    // standard_parallel; defaults to -71
};


PetscErrorCode getFlowLawFromUser(MPI_Comm com, IceType* &ice, PetscInt &flowLawNum);


// this utility prints only when verbosityLevel >= thresh; see iMutil.cc
extern PetscInt verbosityLevel;
PetscErrorCode setVerbosityLevel(PetscInt level);
PetscErrorCode verbosityLevelFromOptions();
PetscErrorCode verbPrintf(const int thresh, MPI_Comm comm,const char format[],...);

#endif /* __iceModelpreamble_hh */
