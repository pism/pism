// Copyright (C) 2007--2010 Nathan Shemonski, Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include <cstring>
#include "../base/pism_const.hh"
#include "../base/nc_util.hh"
#include "forcing.hh"


Data1D::Data1D() {
  vindep = PETSC_NULL;
  vdata = PETSC_NULL;
  interpCode = DATA1D_LINEAR_INTERP;
}


Data1D::~Data1D() {
  if (vindep != PETSC_NULL) {
    VecDestroy(vindep);
  }
  if (vdata != PETSC_NULL) {
    VecDestroy(vdata);
  }
}


//! Read data from NetCDF file (specified by a file name) into a Data1D.
/*! Takes NetCDF file name and names of variables for independent variable
    (e.g. time variable) and dependent variable (e.g. temperature anomaly in
    case of ice core data).  Creates instance of Data1D class.  Reads data on
    processor zero.  Allocates sequential Vec on each processor.  Broadcasts
    processor zero data to all processors.
 */
PetscErrorCode Data1D::readData(MPI_Comm mycom, PetscMPIInt myrank,
                                    const char *myncfilename,
                                    const char *myindepvarname, const char *mydatavarname) {
  PetscErrorCode ierr;
  int stat, ncid = 0;
  if (myrank == 0) {
    stat = nc_open(myncfilename, 0, &ncid);
    if (stat != NC_NOERR) {
      ierr = PetscPrintf(mycom, "Data1D ERROR: Can't open file '%s'.\n", myncfilename); CHKERRQ(ierr);
      PetscEnd();
    }
  }
  MPI_Bcast(&ncid, 1, MPI_INT, 0, mycom);
  ierr = readData(mycom, myrank, ncid, myindepvarname, mydatavarname); CHKERRQ(ierr);
  if (myrank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }
  return 0;
}


//! Read data from NetCDF file (specified by a NetCDF integer id) into a Data1D.
PetscErrorCode Data1D::readData(MPI_Comm mycom, PetscMPIInt myrank,
                                int myncid, 
                                const char *myindepvarname, const char *mydatavarname) {
  PetscErrorCode ierr;
  int indepid, dataid;

  com = mycom;
  rank = myrank;
  strcpy(indepvarname, myindepvarname);
  strcpy(datavarname, mydatavarname);
  if (rank == 0) {
    int stat;
    stat = nc_inq_varid(myncid, indepvarname, &indepid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(myncid, datavarname, &dataid); CHKERRQ(nc_check(stat));
  }
  ierr = getInterpolationCode(myncid, dataid, &interpCode); CHKERRQ(ierr);
  MPI_Bcast(&interpCode, 1, MPI_INT, 0, com);
  ierr = ncVarBcastVec(myncid, indepid, &vindep); CHKERRQ(ierr);  // creates this Vec
  ierr = ncVarBcastVec(myncid, dataid, &vdata); CHKERRQ(ierr);  // creates this Vec
  return 0;
}


PetscErrorCode Data1D::getInterpolationCode(int ncid, int vid, int *code) {
  PetscErrorCode ierr;

  if (rank == 0) {
    int stat;
    char *attr;
    size_t len;

    stat = nc_inq_attlen(ncid, vid, "interpolation", &len);
    if (stat == NC_NOERR) {
      attr = new char[len];
      stat = nc_get_att_text(ncid, vid, "interpolation", attr); CHKERRQ(nc_check(stat));
      attr[len] = '\0';
      if (strcmp(attr, "constant_piecewise_forward") == 0) {
        *code = DATA1D_CONST_PIECE_FWD_INTERP;
      } else if (strcmp(attr, "constant_piecewise_backward") == 0) {
        *code = DATA1D_CONST_PIECE_BCK_INTERP;
      } else if (strcmp(attr, "linear") == 0) {
        *code = DATA1D_LINEAR_INTERP;
      } else {
        ierr = verbPrintf(5, com, 
            "ATTENTION: interpolation '%s' for 1D data %s is unknown; defaulting to linear\n",
            attr,datavarname); CHKERRQ(ierr);
        *code = DATA1D_LINEAR_INTERP;
      }
      delete[] attr;
    } else {
      ierr = verbPrintf(5, com, 
          "ATTENTION: interpolation attribute for 1D data %s is not found; defaulting to linear\n",
          datavarname); CHKERRQ(ierr);
      *code = DATA1D_LINEAR_INTERP;
    }
  } else {
    *code = -1;
  }
  return 0;
}


//! Broadcast sequential Vecs containing 1D climate data to each processor.
PetscErrorCode Data1D::ncVarBcastVec(int ncid, int vid, Vec *vecg) {
  
  PetscErrorCode ierr;
  int stat;
  size_t M_s;
  int M;
  float *f = NULL;

  if (rank == 0) {
    int dimids[NC_MAX_VAR_DIMS];
    int ndims, natts;
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 1) {
      SETERRQ2(1, "number of dimensions = %d for %s; should have ndims=1\n",
               ndims, name);
    }
    stat = nc_inq_dimlen(ncid, dimids[0], &M_s); CHKERRQ(nc_check(stat));
    M = (int)M_s;
    f = new float[M];
    stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
  }
  ierr = MPI_Bcast(&M, 1, MPI_INT, 0, com); CHKERRQ(ierr); // broadcast the length

  if (rank != 0){ // if you're not rank 0, you still need to create the array
    f = new float[M];
  }
  ierr = MPI_Bcast(f, M, MPI_FLOAT, 0, com); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, M, vecg); CHKERRQ(ierr);
  for (int x=0; x<(int)M; x++) {
    ierr = VecSetValue(*vecg, x, f[x], INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(*vecg); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*vecg); CHKERRQ(ierr);

  delete [] f;
  return 0;
}


PetscErrorCode Data1D::getIndexMax(PetscInt *len) {
  PetscErrorCode  ierr = VecGetLocalSize(vdata, len); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode Data1D::getIndexedDataValue(PetscInt index, PetscScalar *value) {
  PetscErrorCode ierr;
  PetscScalar *data;
  PetscInt    len;

  if (index < 0) {
    SETERRQ(1,"index negative");
  }
  ierr = VecGetLocalSize(vdata, &len); CHKERRQ(ierr);
  if (index >= len) {
    SETERRQ(1,"index out of bounds: too large");
  }
  ierr = VecGetArray(vdata, &data); CHKERRQ(ierr);
  *value = data[index];
  ierr = VecRestoreArray(vdata, &data); CHKERRQ(ierr);
  return 0;
}


/*!
Binary search into a 1D array, <tt>double indep[length]</tt>, to return index \c l so that
<tt>indep[l] <= myval <= indep[l+1]</tt>, or so that <tt>l=length-1</tt> and <tt>indep[l] <= myval</tt>,
or so that <tt>l=0</tt> and <tt>myval <= indep[l]</tt>.

In every case, <tt>indep[l]</tt> will be valid, where \c l is the return value, but
when <tt>l=length-1</tt> then <tt>indep[l+1]</tt> is not valid.

It is suggested that calling procedures first check if \c myval is outside the range 
[<tt>indep[0],indep[length-1]</tt>].
 */
PetscInt Data1D::binarySearch(double *indep, const PetscInt length, const double myval) {

  if (myval <= indep[0])           return 0;
  if (myval >= indep[length-1])    return length-1;

  PetscInt l = 0,       // left
           r = length;  // right+1

  while (r > l + 1) {
    PetscInt j = (r + l)/2;
    if (myval < indep[j]) { //     indep[l]  <=  myval  <  indep[j]             <=  indep[r]
      r = j;
    } else {                //     indep[l]  <=            indep[j]  <=  myval  <=  indep[r]
      l = j;
    }
  }

  return l;
}


PetscErrorCode Data1D::getInterpolatedDataValue(PetscScalar myindep, PetscScalar *value) {
  PetscErrorCode ierr;
  PetscScalar *indep, *data;
  PetscInt    index, len;

  // determine index into data
  ierr = VecGetLocalSize(vindep, &len); CHKERRQ(ierr);
  ierr = VecGetArray(vindep, &indep); CHKERRQ(ierr);
  // do a binary search to find where our value of the independent variable fits in;
  //   here we need independent variable to be sorted into increasing order!
  index = binarySearch(indep, len, myindep);
  if (index < 0) {
    SETERRQ(1,"computed index is negative");
  }
  if (index >= len) {
    SETERRQ(2,"computed index exceeds length");
  }
  ierr = VecRestoreArray(vindep, &indep); CHKERRQ(ierr);

  // interpolate
  ierr = VecGetArray(vdata, &data); CHKERRQ(ierr);
  if (myindep == indep[index]) {
    *value = data[index]; // if we have exact data, use it
  } else { // otherwise we need to interpolate
    switch (interpCode) {
      case DATA1D_CONST_PIECE_BCK_INTERP:
        // use the data point we are in front of
        if (index == 0) {
          *value = data[index];
        } else {
          *value = data[index-1];
        }
        break;
      case DATA1D_CONST_PIECE_FWD_INTERP:
        *value = data[index];
        break;
      case DATA1D_LINEAR_INTERP:
        if (index == 0) {
          *value = data[index];
        } else {
          const PetscScalar slope = (myindep - data[index-1]) / (data[index] - data[index-1]);
          *value = data[index-1] + slope * (data[index]- data[index-1]);
        }
        break;
      default:
        SETERRQ1(1, "unknown interpolation method for %s\n", datavarname);
    } // end switch
  }
  ierr = VecRestoreArray(vdata, &data); CHKERRQ(ierr);

  return 0;
}



PISMClimateForcing::PISMClimateForcing() : Data1D() {
  index = 0;
}


//FIXME:  the user of readClimateForcingData() should provide a context (struct) with
//   STRUCT FIELD               EXAMPLE (values for -dTforcing)
//   char     indepname;        "t"
//   char     dataname[30];     "delta_T"
//   char     interpname[30];   "constant_piecewise_forward"
//   char     units[30];        "K"
//   bool     scalebyminus;     true

PetscErrorCode PISMClimateForcing::readClimateForcingData(MPI_Comm mycom, PetscMPIInt myrank,
        int ncid, PetscScalar curr_year, PetscInt datatype) {
  PetscErrorCode ierr;
  
  switch (datatype) {
    case PCF_DELTA_T:
      ierr = readData(mycom,myrank,ncid,"t","delta_T"); CHKERRQ(ierr);
      break;
    case PCF_DELTA_SEA_LEVEL:
      ierr = readData(mycom,myrank,ncid,"t","delta_sea_level"); CHKERRQ(ierr);
      break;
    default:
      SETERRQ(1,"invalid datatype\n");
  } 

  // times are positive ages (= years b.p.) in data; change to negative (years *after* present)
  ierr = VecScale(vindep,-1.0); CHKERRQ(ierr); // puts in increasing order!! FIXME: SILLY!

  ierr = initIndexClimateForcingData(curr_year); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMClimateForcing::initIndexClimateForcingData(PetscScalar curr_year) {
  PetscErrorCode ierr;
  PetscInt       len;
  PetscScalar    *timeinyears;

  ierr = VecGetArray(vindep, &timeinyears); CHKERRQ(ierr);
  ierr = VecGetLocalSize(vindep, &len); CHKERRQ(ierr);

  // do a binary search to find where our current year fits in
  index = binarySearch(timeinyears, len, curr_year);

  ierr = VecRestoreArray(vindep, &timeinyears); CHKERRQ(ierr);

  ierr = verbPrintf(3, com, "index into PISMClimateForcing data %s found using current year: index = %d\n",
          datavarname, index); CHKERRQ(ierr);
  // maybe we are already past our place.
  if (index >= len - 1) {
    ierr = verbPrintf(1, com, 
      "PISMClimateForcing ATTENTION:  At or past end of climate forcing data %s when initializing indexing.\n",
      datavarname); CHKERRQ(ierr);
    index = len - 1; // presumably redundant, but guarantee timeinyears[index] valid
  }
  if (index < 0) { SETERRQ(1,"index negative in PISMClimateForcing::initIndexClimateForcingData()"); }

  return 0;
}


//! Return the interpolated climate data value at time curr_year.
/*!
Consecutive calls to this must use non-decreasing values of curr_year.  An internal index is
advanced at the call, for efficiency.  There is no backward index movement allowed.

If this call is made when PISMClimateForcing::index is out of range, it generates lots of
"PISMClimateForcing ATTENTION: ..." messages at stdout.  Caller beware ...
 */
PetscErrorCode PISMClimateForcing::updateFromClimateForcingData(PetscScalar curr_year, PetscScalar *value) {
  PetscErrorCode ierr;
  PetscScalar    *timeinyears, *data;
  PetscInt       len;

  if (index < 0) { SETERRQ(1, "PISMClimateForcing ERROR: index < 0");  }

  ierr = VecGetArray(vindep, &timeinyears); CHKERRQ(ierr);
  ierr = VecGetLocalSize(vindep, &len); CHKERRQ(ierr);

  // do a binary search to find where our current year fits in
  index = binarySearch(timeinyears, len, curr_year);

  ierr = VecRestoreArray(vindep, &timeinyears); CHKERRQ(ierr);

  if (index < 0) { SETERRQ(2, "PISMClimateForcing ERROR: index < 0");  }

  if (index >= len - 1) {
    ierr = verbPrintf(1, com, 
       "PISMClimateForcing ATTENTION: no more data for climate forcing %s;\n"
       "   curr_year = %f; len = %d; index = %d; returning forcing value 0.0\n", 
       datavarname,curr_year,len,index); CHKERRQ(ierr);
    index = len - 1; // guarantee timeinyears[index] valid for now
    *value = 0.0;
    return 0;
  }

  // if we don't check these cases then we will try "data[index-1]" and get seg fault
  if ( (index == 0) && (    (interpCode == DATA1D_CONST_PIECE_BCK_INTERP)
                         || (interpCode == DATA1D_LINEAR_INTERP)          ) ) {
    ierr = verbPrintf(1, com, 
       "PISMClimateForcing ATTENTION: model year not far enough after beginning of data for interpolation\n"
       "   of climate forcing %s; curr_year = %f; index = %d; returning forcing value 0.0\n",
       datavarname,curr_year,index); CHKERRQ(ierr);
    *value = 0.0;
    return 0;
  }
  
  // so now we should actually look at the data ...
  ierr = VecGetArray(vdata, &data); CHKERRQ(ierr);
  if (curr_year == timeinyears[index]) {
    *value = data[index]; // if we have exact data, use it
  } else { // otherwise we need to interpolate
    switch (interpCode) {
      case DATA1D_CONST_PIECE_BCK_INTERP:
        *value = data[index-1]; // use the data point we are in front of
        break;
      case DATA1D_CONST_PIECE_FWD_INTERP:
        *value = data[index];
        break;
      case DATA1D_LINEAR_INTERP:
        { // scope needed when declaring vars ... C is silly ...
          PetscScalar z0 = data[index-1],
                      z1 = data[index],
                      slope = (z1 - z0) / (timeinyears[index] - timeinyears[index-1]);
          *value = z0 + slope * (curr_year - timeinyears[index-1]);
        }
        break;
      default:
        SETERRQ1(3,"unknown interpolation method in PISMClimateForcing %s\n", datavarname);
    } // end switch
  } // end else
  ierr = VecRestoreArray(vdata, &data); CHKERRQ(ierr);

  return 0;
}

