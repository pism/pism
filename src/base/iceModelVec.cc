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
           "                  long_name = '%s'\n", long_name);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  standard_name = '%s'\n", standard_name);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  units = '%s'\n", units);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  pism_intent = '%s'\n\n", pism_intent);  CHKERRQ(ierr);
  return 0;
}

//! Sets NetCDF attributes of an IceModelVec object.
/*! Call setAttrs("new long name", "new units", "new pism_intent", NULL) if a
  variable does not have a standard name. Similarly, by putting NULL in an
  appropriate spot, it is possible tp leave long_name, units or pism_intent
  unmodified.
 */
PetscErrorCode  IceModelVec::set_attrs(const char my_long_name[], const char my_units[],
				       const char my_pism_intent[], const char my_standard_name[]) {
  if (my_long_name != NULL)
    strcpy(long_name,my_long_name);

  if (my_units != NULL)
    strcpy(units,my_units);

  if (my_pism_intent != NULL)
    strcpy(pism_intent,my_pism_intent);

  if (my_standard_name != NULL) {
    strcpy(standard_name,my_standard_name);
    has_standard_name = PETSC_TRUE;
  }
  return 0;
}

//! Finds the variable by its standard_name attribute, which has to be set using setAttrs
/*!
  Here's how it works:

  1) Check if the current IceModelVec has a standard_name. If it does, go to
  step 2, otherwise go to step 4.

  2) Find the variable with this standard_name. Bail out if two
  variables have the same standard_name, otherwise go to step 3.

  3) If no variable was found, go to step 4, otherwise go to step 5.

  4) Find the variable with the right variable name. Go to step 5.

  5) Broadcast the existence flag and the variable ID.
 */
PetscErrorCode  IceModelVec::find(const int ncid, int *varidp, PetscTruth *exists) {
  PetscInt ierr;
  size_t attlen;
  int stat = 0, found = 0, my_varid = -1, nvars;
  char attribute[TEMPORARY_STRING_LENGTH];

  // Processor 0 does all the job here.
  if (grid->rank == 0) {

    if (has_standard_name) {
      ierr = nc_inq_nvars(ncid, &nvars); CHKERRQ(check_err(stat,__LINE__,__FILE__));

      for (int j = 0; j < nvars; j++) {
	stat = nc_get_att_text(ncid, j, "standard_name", attribute);
	if (stat != NC_NOERR) {
	  continue;
	}

	// text attributes are not always zero-terminated, so we need to add the
	// trailing zero:
	stat = nc_inq_attlen(ncid, j, "standard_name", &attlen);
	CHKERRQ(check_err(stat,__LINE__,__FILE__));
	attribute[attlen] = 0;

	if (strcmp(attribute, standard_name) == 0) {
	  if (!found) {		// if unique
	    found = 1;
	    my_varid = j;
	  } else {    // if not unique
	    fprintf(stderr, "Variables #%d and #%d have the same standard_name ('%s').\n",
		   my_varid, j, attribute);
	    SETERRQ(1,"Inconsistency in the input file: two variables have the same standard_name.");	  
	  }
	}
      }
    } // end of if(has_standard_name)

    if (!found) {
      // look for varname
      stat = nc_inq_varid(ncid, varname, &my_varid);
      if (stat == NC_NOERR)
	found = 1;
    }
  } // end of if(grid->rank == 0)

  // Broadcast the existence flag and the variable ID.
  ierr = MPI_Bcast(&found, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&my_varid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

  if (found) {
    *exists = PETSC_TRUE;
    if (varidp != NULL)
      *varidp = my_varid;
  } else {
    *exists = PETSC_FALSE;
    // *varidp is not modified
  }

  return 0;
}

//! Defines a netcdf variable corresponding to an IceModelVec object. Virtual only.
PetscErrorCode IceModelVec::define_netcdf_variable(int ncid, nc_type nctype, int *varidp) {
  SETERRQ(1, "define_netcdf_variable: virtual only");
}

PetscErrorCode IceModelVec::write_attrs(const int ncid) {
  PetscTruth exists;
  int varid, ierr;
  
  ierr = find(ncid, &varid, &exists); CHKERRQ(ierr);
  if (!exists)
    SETERRQ(1, "Can't write attributes of an undefined variable.");

  if (grid->rank == 0) {
    ierr = nc_put_att_text(ncid, varid,"pism_intent", strlen(pism_intent), pism_intent);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    ierr = nc_put_att_text(ncid, varid,"units", strlen(units), units);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    ierr = nc_put_att_text(ncid, varid,"long_name", strlen(long_name), long_name);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    ierr = nc_put_att_text(ncid, varid,"standard_name", strlen(standard_name), standard_name);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}

//! Calls the appropriate NCTool method to save the IceModelVec into a NetCDF variable.
PetscErrorCode IceModelVec::putVecNC(const int ncid, const int *s, const int *c, int dims, 
                                          void *a_mpi, int a_size) {
  PetscErrorCode ierr;
  NCTool nct(grid);
  PetscTruth exists;
  int varid_nc;

  ierr = find(ncid, &varid_nc, &exists); CHKERRQ(ierr);
  if (!exists)
    SETERRQ(1, "Variable does not exist");

  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.put_local_var(ncid, varid_nc, da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.put_global_var(ncid, varid_nc, da, v,
                          s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  }
  return 0;
}

PetscErrorCode IceModelVec::read(const char filename[], const unsigned int time) {           
  SETERRQ(1, "IceModelVec::read(...) is virtual only.");
  return 0;
}


//! Calls the appropriate NCTool method to read a NetCDF variable into the IceModelVec.
PetscErrorCode IceModelVec::read_from_netcdf(const char filename[], const unsigned int time,
					     int dims, const int Mz) {           
  PetscErrorCode ierr;
  PetscTruth exists;
  void *a_mpi;
  int a_len, max_a_len, ncid, varid;
  NCTool nct(grid);
  int s[] = {time, grid->xs, grid->ys, 0}; // Start local block: t dependent; 
  int c[] = {1, grid->xm, grid->ym, Mz}; // Count local block: t dependent

  // Memory allocation:
  max_a_len = a_len = grid->xm * grid->ym * Mz;
  ierr = MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid->com); CHKERRQ(ierr);
  ierr = PetscMalloc(max_a_len * sizeof(double), &a_mpi); CHKERRQ(ierr);

  if (grid->rank == 0) {
    ierr = nc_open(filename, NC_NOWRITE, &ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  ierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

  ierr = find(ncid, &varid, &exists); CHKERRQ(ierr);
  if (!exists)
    SETERRQ2(1, "Can't find '%s' in '%s'.", varname, filename);

  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.get_local_var_id(ncid, varid, da, v, g,
                         s, c, dims, a_mpi, max_a_len); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.get_global_var_id(ncid, varid, da, v,
                          s, c, dims, a_mpi, max_a_len); CHKERRQ(ierr);  
  }

  if (grid->rank == 0) {
    ierr = nc_close(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  ierr = PetscFree(a_mpi);
  return 0;
}

//! Writes an IceModelVec to a NetCDF file.
/*!
  1) Get the last time.
  2) Find the variable in the file. Call define_variable if not found.
  3) Call put_global_var or put_local_var.
 */
PetscErrorCode IceModelVec::write_to_netcdf(const int ncid, int dims, nc_type nctype,
					    const int Mz, void *a_mpi, int a_size) {
  PetscErrorCode ierr;
  PetscTruth exists;
  NCTool nct(grid);
  int t_id, varid;
  int t;

  // get the last time (index, not the value):
  if (grid->rank == 0) {
    size_t t_len;
    ierr = nc_inq_dimid(ncid, "t", &t_id); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    ierr = nc_inq_dimlen(ncid, t_id, &t_len); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    t = (int)t_len;
  }
  ierr = MPI_Bcast(&t, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  int s[] = {t - 1, grid->xs, grid->ys, 0}; // Start local block: t dependent; 
  int c[] = {1, grid->xm, grid->ym, Mz}; // Count local block: t dependent

  // find the variable
  ierr = find(ncid, &varid, &exists); CHKERRQ(ierr);
  if (!exists) {
    ierr = define_netcdf_variable(ncid, nctype, &varid); CHKERRQ(ierr);
  }

  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.put_local_var(ncid, varid, da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.put_global_var(ncid, varid, da, v,
                          s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  }
  return 0;
  
}


//! Calls the appropriate NCTool method to regrid a NetCDF variable from some file into the IceModelVec.
PetscErrorCode  IceModelVec::regridVecNC(int dim_flag, LocalInterpCtx &lic)  {
  PetscErrorCode ierr;
  NCTool nct(grid);
  // FIXME: a flag for whether the IceModelVec is really integer-valued should be checked; if so 
  // then regrid_local|global_var() should be called w last arg "true", after setting the MaskInterp
  // for NCTool
  if (localp) {
    Vec g;
    ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
    ierr = nct.regrid_local_var(varname, dim_flag, lic, da, v, g, false); CHKERRQ(ierr);
    ierr = VecDestroy(g); CHKERRQ(ierr);
  } else {
    ierr = nct.regrid_global_var(varname, dim_flag, lic, da, v, false); CHKERRQ(ierr);
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
