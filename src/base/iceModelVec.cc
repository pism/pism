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
#include "pism_const.hh"
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
}


IceModelVec::~IceModelVec() {
}


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


PetscErrorCode  IceModelVec::printInfo(const PetscInt verbosity) {
  PetscErrorCode ierr;
  
  if (grid == PETSC_NULL) {
    SETERRQ1(1,"ERROR: cannot print info for IceModelVec with varname='%s'\n"
               "   because grid=PETSC_NULL.  ENDING.\n\n",varname);
  }

  ierr = verbPrintf(verbosity,grid->com,
         "\nprinting info for IceModelVec with varname='%s':\n",
         varname); CHKERRQ(ierr);
  if (da == PETSC_NULL) {
    ierr = verbPrintf(verbosity,grid->com,
          "  WARNING:  da == PETSC_NULL for IceModelVec with varname='%s'!\n",
          varname); CHKERRQ(ierr);
  }
  if (v == PETSC_NULL) {
    ierr = verbPrintf(verbosity,grid->com,
          "  WARNING:  v == PETSC_NULL for IceModelVec with varname='%s'!\n",
          varname); CHKERRQ(ierr);
  }
  if (array == PETSC_NULL) {
    ierr = verbPrintf(verbosity,grid->com,
          "  WARNING:  array == PETSC_NULL for IceModelVec with varname='%s'!\n",
          varname); CHKERRQ(ierr);
  }
  
  ierr = verbPrintf(verbosity,grid->com,
           "  boolean flags:  localp = %d,  IOwnDA = %d,  has_standard_name = %d\n",
           (int)localp, (int)IOwnDA, has_standard_name);  CHKERRQ(ierr);

  ierr = verbPrintf(verbosity,grid->com,
           "  NetCDF info:    varid_nc = %d\n", varid_nc);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  long_name = '%s'\n", long_name);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  standard_name = '%s'\n", standard_name);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  units = '%s'\n", units);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  pism_intent = '%s'\n\n", pism_intent);  CHKERRQ(ierr);
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
PetscErrorCode IceModelVec::getVecNC(const int ncid, const int *s, const int *c, 
                             int dims, void *a_mpi, int a_size) {           
  PetscErrorCode ierr;
  NCTool nct;
  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.get_local_var(grid, ncid, varname, da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.get_global_var(grid, ncid, varname, da, v,
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
    ierr = nct.put_local_var(grid, ncid, varid_nc, da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.put_global_var(grid, ncid, varid_nc, da, v,
                          s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  }
  return 0;
}


//! Calls the appropriate NCTool method to regrid a NetCDF variable from some file into the IceModelVec.
PetscErrorCode  IceModelVec::regridVecNC(int dim_flag, LocalInterpCtx &lic)  {
  PetscErrorCode ierr;
  NCTool nct;
  // FIXME: a flag for whether the IceModelVec is really integer-valued should be checked; if so 
  // then regrid_local|global_var() should be called w last arg "true", after setting the MaskInterp
  // for NCTool
  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.regrid_local_var(varname, dim_flag, lic, *grid, da, v, g, false); CHKERRQ(ierr);
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.regrid_global_var(varname, dim_flag, lic, *grid, da, v, false); CHKERRQ(ierr);
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


/********* IceModelVec3 and IceModelVec3Bedrock: SEE SEPARATE FILE  iceModelVec3.cc    **********/

/********* IceModelVec2 and IceModelVec2Box: SEE SEPARATE FILE  iceModelVec2.cc    **********/
