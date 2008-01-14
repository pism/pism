// Copyright (C) 2008 Ed Bueler
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

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"

#include "iceModelVec.hh"


IceModelVec::IceModelVec() {

  v = PETSC_NULL;
  da = PETSC_NULL;
  grid = PETSC_NULL;
  array = PETSC_NULL;
  localp = true;
  IOwnDA = true;

  strcpy(varname,"*****UNKNOWN**** varname (NetCDF or not)");

  strcpy(long_name,"UNKNOWN NetCDF long_name");
  strcpy(units,"UNKNOWN NetCDF units");
  strcpy(pism_intent,"UNKNOWN NetCDF pism_intent");

  has_standard_name = PETSC_FALSE;
  strcpy(standard_name,"UNKNOWN NetCDF CF 1.0 standard_name");
  
  varid_nc = -9999;
};


IceModelVec::~IceModelVec() {
};


PetscErrorCode  IceModelVec::create(IceGrid &mygrid, const char my_varname[], bool local) {
  SETERRQ(1,"VIRTUAL ONLY: not implemented");
  return 0;
}


PetscErrorCode  IceModelVec::destroy() {
  PetscErrorCode ierr;
  if (v != PETSC_NULL) {
    ierr = VecDestroy(v); CHKERRQ(ierr);
    v = PETSC_NULL;
  }
  if ((IOwnDA) && (da != PETSC_NULL)) {
    ierr = DADestroy(da); CHKERRQ(ierr);
    da = PETSC_NULL;
  }
  return 0;
}


PetscErrorCode  IceModelVec::setVaridNC(const int my_varid) {
  varid_nc = my_varid;
  return 0;
}


PetscErrorCode  IceModelVec::setAttrsNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[]) {
  varid_nc = my_varid;
  strcpy(long_name,my_long_name);
  strcpy(units,my_units);
  strcpy(pism_intent,my_pism_intent);
  has_standard_name = PETSC_FALSE;
  return 0;
}


PetscErrorCode  IceModelVec::setAttrsCFstandardNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[],
             const char my_standard_name[]) {
  varid_nc = my_varid;
  strcpy(long_name,my_long_name);
  strcpy(units,my_units);
  strcpy(pism_intent,my_pism_intent);
  strcpy(standard_name,my_standard_name);
  has_standard_name = PETSC_TRUE;
  return 0;
}


PetscErrorCode  IceModelVec::writeAttrsNC(const int ncid) {
  SETERRQ(1,"not YET implemented");
  return 0;
}


PetscErrorCode  IceModelVec::findVecNC(const int ncid, PetscTruth *exists) {
  SETERRQ(1,"not YET implemented");
/*
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // on processor zero, use varname string to check if it is present
  if (grid->rank == 0) {
    int stat = nc_inq_varid(ncid, varname_nc, &varid_nc); 
    exists = (stat == NC_NOERR);
  }
  // broadcast the existence flag
  ierr = MPI_Bcast(&exists, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  // if it exists, read it
  if (exists == PETSC_TRUE) {
    ierr = ncVarToDAVec(ncid, varid_nc, *da, v, g2, vzero); CHKERRQ(ierr);
  }
*/
  return 0;
}


//! Calls the appropriate NCTool method to read a NetCDF variable into the IceModelVec.
PetscErrorCode IceModelVec::getVecNC(const int ncid, const int *s, const int *c, int dims, 
                                     void *a_mpi, int a_size) {           
  PetscErrorCode ierr;
  NCTool nct;
  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.get_local_var(grid, ncid, varname, NC_FLOAT, da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.get_global_var(grid, ncid, varname, NC_FLOAT, da, v,
                          s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  }
  return 0;
}


//! Calls the appropriate NCTool method to save the IceModelVec into a NetCDF variable.
PetscErrorCode IceModelVec::putVecNC(const int ncid, const int *s, const int *c, int dims, 
                                          void *a_mpi, int a_size) {
  PetscErrorCode ierr;
  NCTool nct;
  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.put_local_var(grid, ncid, varid_nc, NC_FLOAT, da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.put_global_var(grid, ncid, varid_nc, NC_FLOAT, da, v,
                          s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  }
  return 0;
}


//! Calls the appropriate NCTool method to regrid a NetCDF variable from some file into the IceModelVec.
PetscErrorCode  IceModelVec::regridVecNC(const char *vars, char c, int dim_flag, LocalInterpCtx &lic)  {
  PetscErrorCode ierr;
  NCTool nct;
  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.regrid_local_var(vars, c, varname, dim_flag, lic, *grid, da, v, g); CHKERRQ(ierr);
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.regrid_global_var(vars, c, varname, dim_flag, lic, *grid, da, v); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode  IceModelVec::checkAllocated() {
  if (v == PETSC_NULL) {
    SETERRQ1(1,"IceModelVec ERROR: IceModelVec with varname='%s' NOT allocated\n",
             varname);
  }
  return 0;
}


PetscErrorCode  IceModelVec::checkHaveArray() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  if (array == PETSC_NULL) {
    SETERRQ1(1,"array for IceModelVec with varname='%s' not available\n"
               "  (REMEMBER TO RUN needAccessToVals() before access and doneAccessToVals() after access)\n",
               varname);
  }
  return 0;
}


/*
PetscErrorCode  IceModelVec::checkSelfOwnsIt(const PetscInt i, const PetscInt j) {
  if (allocated == PETSC_FALSE) {
    SETERRQ3(1,"IceModelVec ERROR: (i,j)=(%d,%d) not in ownership range of processor %d\n",
             i,j,grid->rank);
  }
  return 0;
}
*/


PetscErrorCode  IceModelVec::needAccessToVals() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DAVecGetArray(da, v, &array); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::doneAccessToVals() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da, v, &array); CHKERRQ(ierr);
  array = PETSC_NULL;
  return 0;
}


PetscErrorCode  IceModelVec::beginGhostComm() {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has varname='%s')\n",
               varname);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(da, v, INSERT_VALUES, v);  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::endGhostComm() {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has varname='%s')\n",
               varname);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(da, v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::setToConstant(const PetscScalar c) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = VecSet(v,c); CHKERRQ(ierr);
  return 0;
}



/********* IceModelVec2 **********/

IceModelVec2::IceModelVec2() : IceModelVec() {
};


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with varname='%s' already allocated\n",my_varname);
  }
  PetscErrorCode ierr = create(my_grid, my_varname, local, DA_STENCIL_STAR); CHKERRQ(ierr);
  return 0;
}
  

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_varname[], bool local,
                                     DAStencilType my_sten) {

  grid = &my_grid;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate2d(my_grid.com, DA_XYPERIODIC, my_sten, N, M, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
  IOwnDA = true;
  
  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(varname,my_varname);
  return 0;
}


// Return value of ice scalar quantity stored in an IceModelVec2.
PetscScalar     IceModelVec2::getVal(const PetscInt i, const PetscInt j) {
  checkHaveArray();
  PetscScalar **arr = (PetscScalar**) array;
  return arr[i][j];
}


// Return values on planar star stencil of ice scalar quantity stored in an IceModelVec2.
PetscErrorCode   IceModelVec2::getPlaneStar(const PetscInt i, const PetscInt j, planeStar *star) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec2 ERROR: IceModelVec2 with varname='%s' is GLOBAL and cannot do getPlaneStar()\n",
             varname);
  }
  
  PetscScalar **arr = (PetscScalar**) array;

  star->ij = arr[i][j];
  star->ip1 = arr[i+1][j];
  star->im1 = arr[i-1][j];
  star->jp1 = arr[i][j+1];
  star->jm1 = arr[i][j-1];
  return 0;
}


PetscScalar**   IceModelVec2::arrayGet() {
  checkHaveArray();
  return (PetscScalar**) array;
}



/********* IceModelVec2Box **********/

IceModelVec2Box::IceModelVec2Box() : IceModelVec2() {
};


PetscErrorCode  IceModelVec2Box::create(IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2Box with varname='%s' already allocated\n",my_varname);
  }
  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_varname, local, DA_STENCIL_BOX); CHKERRQ(ierr);
  return 0;
}


// Return values on planar BOX stencil of ice scalar quantity stored in an IceModelVec2Box.
PetscErrorCode   IceModelVec2Box::getPlaneBox(const PetscInt i, const PetscInt j, planeBox *box) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec2Box ERROR: IceModelVec2Box with varname='%s' is GLOBAL and cannot do getPlaneBox()\n",
             varname);
  }

  PetscScalar **arr = (PetscScalar**) array;

  box->ij     = arr[i][j];
  box->ip1    = arr[i+1][j];
  box->im1    = arr[i-1][j];
  box->jp1    = arr[i][j+1];
  box->jm1    = arr[i][j-1];
  box->ip1jp1 = arr[i+1][j+1];
  box->im1jp1 = arr[i-1][j+1];
  box->ip1jm1 = arr[i+1][j-1];
  box->im1jm1 = arr[i-1][j-1];

  return 0;
}


/********* IceModelVec3Bedrock **********/

IceModelVec3Bedrock::IceModelVec3Bedrock() : IceModelVec() {
};


//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3Bedrock::create(IceGrid &my_grid, const char my_varname[], bool local) {

  strcpy(varname,my_varname);

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3Bedrock with varname='%s' already allocated\n",varname);
  }
  if (local) {
    SETERRQ1(2,"IceModelVec3Bedrock must be GLOBAL (varname='%s')\n",varname);
  }

  grid = &my_grid;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(my_grid.com, DA_YZPERIODIC, DA_STENCIL_STAR, my_grid.p->Mbz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
  IOwnDA = true;

  ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);

  localp = false;
  return 0;
}


//! Set values of bedrock scalar quantity at internal levels determined by IceGrid.
/*!
Array \c valsIN must be an allocated array of \c (grid->p)->Mbz \c PetscScalar s.
 */
PetscErrorCode  IceModelVec3Bedrock::setInternalColumn(const PetscInt i, const PetscInt j, 
                                                       PetscScalar *valsIN) {
  
  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k = 0; k < (grid->p)->Mbz; k++) {
    arr[i][j][k] = valsIN[k];
  }
  return 0;
}


//! Set values of bedrock scalar quantity: set all values in a column to the same value.
PetscErrorCode  IceModelVec3Bedrock::setToConstantColumn(
                        const PetscInt i, const PetscInt j, const PetscScalar c) {

  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k = 0; k < (grid->p)->Mbz; k++) {
    arr[i][j][k] = c;
  }
  return 0;
}


//! Return values of bedrock scalar quantity at internal levels determined by IceGrid.
/*!
Return array \c valsOUT is an allocated array of \c (grid->p)->Mbz \c PetscScalar s.
 */
PetscErrorCode  IceModelVec3Bedrock::getInternalColumn(const PetscInt i, const PetscInt j, 
                                                       PetscScalar **valsOUT) {
  
  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  *valsOUT = arr[i][j];
  return 0;
}


/********* IceModelVec3:    SEE SEPARATE FILE  iceModelVec3.cc    **********/

