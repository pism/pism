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
/*!
This class is a general facility for reading one-dimensional data from a NetCDF file
and putting a copy of it on each processor, and for accessing it either by integer index
or by giving the value of the independent variable and interpolating.

There are two vectors of numbers, \c vindep is the independent variable (frequently
the time axis) and \c vdata is the dependent variable.  Thus if t=\c vindep and z=\c vdata
then z=z(t).

The independent variable \c vindep is assumed to be in increasing order when the method
getInterpolatedDataValue() gets called, because that method does binary search.
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
  PetscMPIInt  rank;
 
  PetscInt     binarySearch(double *indep, const PetscInt length, const double myval);
 
private:
  PetscErrorCode getInterpolationCode(int ncid, int vid, int *code);
  PetscErrorCode ncVarBcastVec(int ncid, int vid, Vec *vecg);
};


// codes for datatype in call to PISMClimateForcing::readClimateForcingData()
#define PCF_DELTA_T          0
#define PCF_DELTA_SEA_LEVEL  1

//! Class for reading and storing climate data, like ice or seabed core data, from NetCDF file, and providing it to a PISMClimateCoupler.
/*! 
This derived class of Data1D maintains an additional index into the time series data
based on a current year.  The method updateFromClimateForcingData() returns the value of
the climate variable.  Usually this value is a change from present.

The vindep member of Data1D is interpreted as a time in years.

FIXME:  THIS COMMENT IS ACCURATE BUT THE CHOICES ARE SILLY:
The semantics of the independent variable vindep are that vindep must be increasing for the calls
to the Data1D methods to work.  But it is assumed that the read NetCDF data has positive ages
(i.e. years before present).  So vindep is the t variable of the NetCDF file, *but multiplied
by -1.0*.  Thus vindep is a bunch of negative real numbers in increasing order, increasing
as real numbers.  For example, in the grip case,
   indep[0] = -249900.0, indep[1] = -249850.0, indep[2] = -249800.0, ...
        indep[N-3] = -100.0, indep[N-2] = -50, indep[N-1] = -0.0
while in the NetCDF file we see
   t = 249900, 249850, 249800, ..., 100, 50, 0
 */
class PISMClimateForcing : public Data1D {
public:
  PISMClimateForcing();
  
  PetscErrorCode readClimateForcingData(MPI_Comm mycom, PetscMPIInt myrank,
                                     int ncid, PetscScalar curr_year, PetscInt datatype);
  PetscErrorCode updateFromClimateForcingData(PetscScalar curr_year, PetscScalar *value);

protected:
  PetscInt     index;

private:
  PetscErrorCode initIndexClimateForcingData(PetscScalar curr_year);
};

#endif /* __forcing_hh */

