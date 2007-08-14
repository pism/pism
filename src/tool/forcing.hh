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

#ifndef __forcing_hh
#define __forcing_hh

// codes for interpolation
#define DATA1D_CONST_PIECE_FWD_INTERP  0
#define DATA1D_CONST_PIECE_BCK_INTERP  1
#define DATA1D_LINEAR_INTERP           2

//! Class for reading and storing one-dimensional data on each processor.
/*! This class is a general facility for reading one-dimensional data from a
    NetCDF file and putting a copy of it on each processor, and for accessing it
    either by integer index or by giving the value of the independent variable and
    interpolating.
 */
class Data1D {
public:
  Data1D();
  ~Data1D();
  // first readData() opens and closes NetCDF file
  PetscErrorCode readData(MPI_Comm mycom, PetscMPIInt myrank,
                          const char *myncfilename,
                          const char *myindepvarname, const char *mydatavarname);
  // second readData() assumes NetCDF file is open and does not close it  
  PetscErrorCode readData(MPI_Comm mycom, PetscMPIInt myrank,
                          int myncid,
                          const char *myindepvarname, const char *mydatavarname);
  PetscErrorCode getIndexMax(PetscInt *len);  
  // in getIndexedDataValue(), use  index = 0,1,...,len-1   if len is from getIndexMax()
  PetscErrorCode getIndexedDataValue(PetscInt index, PetscScalar *value);
  PetscErrorCode getInterpolatedDataValue(PetscScalar myindep, PetscScalar *value);

protected:
  Vec          vindep, vdata;
  char         indepvarname[PETSC_MAX_PATH_LEN], datavarname[PETSC_MAX_PATH_LEN];
  PetscInt     interpCode;
  MPI_Comm     com;
  
private:
  PetscErrorCode getInterpolationCode(int ncid, int vid, int *code);
  PetscErrorCode ncVarBcastVec(int ncid, int vid, Vec *vecg);
  PetscMPIInt    rank;
  PetscTruth     vecsAllocated;
};


// codes for datatype in call to IceSheetForcing::readStandardClimateData()
#define ISF_DELTA_T          0
#define ISF_DELTA_SEA_LEVEL  1

class IceSheetForcing : public Data1D {
public:
  IceSheetForcing();
  ~IceSheetForcing();
  
  PetscErrorCode readStandardIceCoreClimateData(MPI_Comm mycom, PetscMPIInt myrank,
                             int ncid, PetscScalar curr_year, PetscInt datatype);
  PetscErrorCode updateFromStandardIceCoreData(PetscScalar curr_year, PetscScalar *change);

protected:
  PetscErrorCode initStandardIceCoreIndex(PetscScalar curr_year);
  PetscInt     index;
  Vec          vtimeinyears;
  PetscTruth   forcingActive;
};

#endif /* __forcing_hh */

