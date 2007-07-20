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

#include <petscda.h>
#include <cstring>
#include <netcdf.h>
#include "forcing.hh"
// next is only needed for verbPrintF()
#include "iceModel.hh"


PetscErrorCode nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}


IceSheetForcing::IceSheetForcing() {
  vecsAllocated = PETSC_FALSE;
  forcingActive = PETSC_FALSE;
  index = 0;
  interpCode = ISF_LINEAR_INTERP;
}


IceSheetForcing::~IceSheetForcing() {
  if (vecsAllocated == PETSC_TRUE) {
    VecDestroy(vtimeinyears);
    VecDestroy(vdata);
    vecsAllocated = PETSC_FALSE;
  }
  forcingActive = PETSC_FALSE;
}


PetscErrorCode IceSheetForcing::readDataAndAlloc(
                           MPI_Comm mycom, PetscMPIInt myrank, PetscInt mydatatype,
                           int ncid, PetscScalar curr_year) {
  PetscErrorCode ierr;
  int stat, timeid, dataid;

  com = mycom;
  rank = myrank;
  datatype = mydatatype;
  switch (datatype) {
    case ISF_DELTA_T:
      strcpy(datavarname, "delta_T");
      break;
    case ISF_DELTA_SEA_LEVEL:
      strcpy(datavarname, "delta_sea_level");
      break;
    default:
      verbPrintf(1,com, "invalid datatype when setting object IceSheetForcing\n");
      PetscEnd();
  } 
  if (rank == 0) {
    stat = nc_inq_varid(ncid, "t", &timeid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, datavarname, &dataid); CHKERRQ(nc_check(stat));
    ierr = getInterpolationCode(ncid, dataid, &interpCode); CHKERRQ(ierr);
  }
  MPI_Bcast(&interpCode, 1, MPI_INT, 0, com);
  ierr = ncVarBcastVec(ncid, timeid, &vtimeinyears); CHKERRQ(ierr);  // creates this Vec
  ierr = ncVarBcastVec(ncid, dataid, &vdata); CHKERRQ(ierr);  // creates this Vec
  ierr = VecScale(vtimeinyears,-1.0); CHKERRQ(ierr);  // vtime values now negative; years *after* present
  vecsAllocated = PETSC_TRUE;
  ierr = initIndex(curr_year); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceSheetForcing::ncVarBcastVec(int ncid, int vid, Vec *vecg) {
  // spread (broadcast) sequential Vecs containing ice or sea bed core-derived 
  // climate data to each processor
  
  PetscErrorCode ierr;
  int stat;
  size_t M;
  float *f = NULL;

  if (rank == 0) {
    int dimids[NC_MAX_VAR_DIMS];
    int ndims, natts;
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 1) {
      SETERRQ2(1, "ncVarBcastDaVec: number of dimensions = %d for %s; should have ndims=1\n",
               ndims, name);
    }
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));
    f = new float[M];
    stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
  }
  ierr = MPI_Bcast(&M, 1, MPI_INT, 0, com); CHKERRQ(ierr); // broadcast the length
 
  // if you're not rank 0, you still need to create the array
  if (rank != 0){
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


PetscErrorCode IceSheetForcing::initIndex(PetscScalar curr_year) {
  PetscErrorCode ierr;
  PetscInt    len;
  PetscScalar *timeinyears;

  ierr = VecGetArray(vtimeinyears, &timeinyears); CHKERRQ(ierr);
  ierr = VecGetLocalSize(vtimeinyears, &len); CHKERRQ(ierr);

  int r, l=0;
  r = len;
  // do a binary search to find where our year fits in.
  while (r > l + 1) {
    PetscInt j = (r + l)/2;
    if(curr_year < timeinyears[j]) {
      r = j;
    } else {
      l = j;
    }
  }    
  index = l;
  // maybe we are already past our place.
  if (l >= len) {
    forcingActive = PETSC_FALSE;
    ierr = verbPrintf(1, com, 
             "ATTENTION: past end of climate forcing data %s.  Using last value.\n",datavarname);
             CHKERRQ(ierr);
  } else {
    forcingActive = PETSC_TRUE;
  }

  ierr = VecRestoreArray(vtimeinyears, &timeinyears); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode IceSheetForcing::getInterpolationCode(int ncid, int vid, int *code) {
  PetscErrorCode ierr;
  int stat;
  char attr[NC_MAX_NAME+1];
  size_t len;

  stat = nc_get_att_text(ncid, vid, "interpolation", attr); CHKERRQ(nc_check(stat));
  stat = nc_inq_attlen(ncid, vid, "interpolation", &len); CHKERRQ(nc_check(stat));
  attr[len] = '\0';
  if (strcmp(attr, "constant_piecewise_forward") == 0) {
    *code = ISF_CONST_PIECE_FWD_INTERP;
  } else if (strcmp(attr, "constant_piecewise_backward") == 0) {
    *code = ISF_CONST_PIECE_BCK_INTERP;
  } else if (strcmp(attr, "linear") == 0) {
    *code = ISF_LINEAR_INTERP;
  } else {
    ierr = verbPrintf(1, com, 
            "ATTENTION: Interpolation '%s' for climate forcing %s is unknown, defaulting to linear.\n",
            attr,datavarname); CHKERRQ(ierr);
    *code = ISF_LINEAR_INTERP;
  }
  return 0;
}


PetscErrorCode IceSheetForcing::updateFromData(PetscScalar curr_year, PetscScalar *change) {
  PetscErrorCode ierr;
  PetscScalar *timeinyears, *data;
  PetscInt    len;

  ierr = VecGetArray(vtimeinyears, &timeinyears); CHKERRQ(ierr);
  ierr = VecGetLocalSize(vtimeinyears, &len); CHKERRQ(ierr);
  ierr = VecGetArray(vdata, &data); CHKERRQ(ierr);

  // if there was a large time step, it is possible
  // that we skip over multiple entries
  while (index < len && curr_year > data[index]) {
    index++;
  }
  if (index >= len) {
    ierr = verbPrintf(1, com, 
             "ATTENTION: no more data for climate forcing %s.\n", datavarname); CHKERRQ(ierr);
    forcingActive = PETSC_FALSE;
  } else {
    if (curr_year == timeinyears[index]) {
      *change = data[index]; // if we have exact data, use it
    } else { // otherwise we need to interpolate
      PetscScalar y0, y1;
      switch (interpCode) {
        case ISF_CONST_PIECE_BCK_INTERP:
          // use the data point we are infront of
          if (index == 0) {
            *change = 0.0;
            ierr = verbPrintf(1, com, 
                     "ATTENTION: model year precedes beginning of data for climate forcing %s;"
                     " setting change=0\n", datavarname); CHKERRQ(ierr);
          } else {
            *change = data[index-1];
          }
          break;
        case ISF_CONST_PIECE_FWD_INTERP:
          *change = data[index];
          break;
        case ISF_LINEAR_INTERP:
          if (index == 0) {
            *change = 0.0;
            ierr = verbPrintf(1, com, 
                     "ATTENTION: model year precedes beginning of data for climate forcing %s;"
                     " setting change=0\n", datavarname); CHKERRQ(ierr);
          } else {
            y0 = data[index-1];
            y1 = data[index];
            *change = y0 + ((y1-y0) / (timeinyears[index] - timeinyears[index-1]))
                                  * (curr_year - timeinyears[index-1]);
          }
          break;
        default:
          SETERRQ1(1, "Unknown interpolation method for climate forcing %s\n", datavarname);
      } // end switch
    }
  }
  
  ierr = VecRestoreArray(vtimeinyears, &timeinyears); CHKERRQ(ierr);
  ierr = VecRestoreArray(vdata, &data); CHKERRQ(ierr);
  return 0;
}

