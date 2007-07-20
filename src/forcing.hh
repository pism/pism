// Copyright (C) 2007 Nathan Shemonski and Ed Bueler
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

#include <petscvec.h>

/* 
[[These utilities are intended to grow into a more general climate forcing facility.]]

For now, when using pgrn, for each of the -dTforcing and -dSLforcing options, a NetCDF file is read.
1.  If the -dTforcing option is set then a variable delta_T is expected.  That (one dimensional)
    data will be used to offset the surface temperature everywhere.
2.  If the -dSLforcing option is set the a variable delta_sea_level is expected.  That
    1D data will be used to move the bed elevations up and down accordingly.
*/

// codes for datatype:
#define ISF_DELTA_T          0
#define ISF_DELTA_SEA_LEVEL  1

// codes for interpolation
#define ISF_CONST_PIECE_FWD_INTERP  0
#define ISF_CONST_PIECE_BCK_INTERP  1
#define ISF_LINEAR_INTERP           2

#ifndef __forcing_hh
#define __forcing_hh

class IceSheetForcing {
public:
  IceSheetForcing();
  ~IceSheetForcing();
  PetscErrorCode readDataAndAlloc(MPI_Comm mycom, PetscMPIInt myrank, PetscInt mydatatype, 
                                  int ncid, PetscScalar curr_year);
  PetscErrorCode updateFromData(PetscScalar curr_year, PetscScalar *change);
  PetscErrorCode printCurrData(PetscScalar curr_year);

protected:
  PetscErrorCode initIndex(PetscScalar curr_year);
  PetscErrorCode getInterpolationCode(int ncid, int vid, int *code);
  PetscErrorCode ncVarBcastVec(int ncid, int vid, Vec *vecg);

private:
  MPI_Comm     com;
  PetscMPIInt  rank;
  PetscInt     datatype, index, interpCode;
  char         datavarname[PETSC_MAX_PATH_LEN];
  Vec          vtimeinyears, vdata;
  PetscTruth   vecsAllocated, forcingActive;
};


#endif /* __forcing_hh */
