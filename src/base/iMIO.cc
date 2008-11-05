// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"


bool IceModel::hasSuffix(const char* fname, const char *suffix) const {
  int flen = strlen(fname);
  int slen = strlen(suffix);
  if (strcmp(fname + flen - slen, suffix) == 0) {
    return true;
  } else {
    return false;
  }
}


PetscErrorCode  IceModel::setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID) {
  PetscErrorCode ierr;

  // read options about year of start, year of end, number of run years;
  // note grid.year has already been set from input file
  PetscScalar usrStartYear, usrEndYear, usrRunYears;
  PetscTruth ysSet, yeSet, ySet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ys", &usrStartYear, &ysSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ye", &usrEndYear, &yeSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-y", &usrRunYears, &ySet); CHKERRQ(ierr);
  if (ysSet == PETSC_TRUE) {
    // user option overwrites data
    ierr = setStartYear(usrStartYear); CHKERRQ(ierr);
    grid.year = usrStartYear;
  } else if (grid_p_year_VALID == PETSC_TRUE) {
    ierr = setStartYear(grid.year); CHKERRQ(ierr);
  } // else do nothing; defaults are set
  if (yeSet == PETSC_TRUE) {
    if (usrEndYear < startYear) {
      SETERRQ(1,
        "ERROR: -ye value less than -ys value (or input file year or default).\n"
        "PISM cannot run backward in time");
    }
    if (ySet == PETSC_TRUE) {
      ierr = verbPrintf(1,grid.com,"WARNING: -y option ignored.  -ye used instead.\n"); CHKERRQ(ierr);
    }
    ierr = setEndYear(usrEndYear); CHKERRQ(ierr);
  } else if (ySet == PETSC_TRUE) {
    ierr = setEndYear(usrRunYears + startYear); CHKERRQ(ierr);
  } else {
    ierr = setEndYear(run_year_default + startYear); CHKERRQ(ierr);
  }
  
  yearsStartRunEndDetermined = PETSC_TRUE;
  return 0;
}

  
PetscErrorCode  IceModel::writeFiles(const char* defaultbasename) {
  PetscErrorCode ierr = writeFiles(defaultbasename,PETSC_FALSE); CHKERRQ(ierr);
  return 0;
}

  
//! Save model state in NetCDF format (and save variables in Matlab format if desired).
/*! 
Optionally allows saving of full velocity field.

Calls dumpToFile_netCDF(), dumpToFile_diagnostic_netCDF(), and writeMatlabVars() to do 
the actual work.
 */
PetscErrorCode  IceModel::writeFiles(const char* defaultbasename,
                                     const PetscTruth forceFullDiagnostics) {
  PetscErrorCode ierr;

  char b[PETSC_MAX_PATH_LEN];
  char fmt[PETSC_MAX_PATH_LEN] = "n";
  char ncf[PETSC_MAX_PATH_LEN]; // netCDF format

  if (doPDD == PETSC_TRUE) { // want to save snow accumulation map, not net accumulation
    ierr = putBackSnowAccumPDD(); CHKERRQ(ierr);
  }
  
  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  // Use the defaults passed from the driver if not specified on command line.
  // We should leave space for a suffix and null byte
  strncpy(b, defaultbasename, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-of", fmt, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);

  if (strchr(fmt, 'p') != NULL) {
    strcat(b,"_pb");  // will write basename_pb.nc
    strcat(fmt,"n");
    ierr = verbPrintf(1, grid.com, 
       "\nWARNING: .pb format no longer supported; writing to NetCDF file %s.nc\n",b);
       CHKERRQ(ierr);
  }

  if (strchr(fmt, 'n') != NULL) {
    strcpy(ncf, b);
    strcat(ncf, ".nc");
    PetscTruth userWantsFull;
    ierr = PetscOptionsHasName(PETSC_NULL, "-f3d", &userWantsFull); CHKERRQ(ierr);
    if ((forceFullDiagnostics == PETSC_TRUE) || (userWantsFull == PETSC_TRUE)) {
      ierr = verbPrintf(2, grid.com, 
            "Writing model state, with full 3D velocities, to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_diagnostic_netCDF(ncf); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com, "Writing model state to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_netCDF(ncf); CHKERRQ(ierr);
    }
  }

  if (strchr(fmt, 'm') != NULL) {
    ierr = verbPrintf(1, grid.com, 
       "\nWARNING: .m format no longer supported with '-o'; use '-mato' and '-matv'\n");
       CHKERRQ(ierr);
  }

  // write out individual variables out to Matlab file
  char       matf[PETSC_MAX_PATH_LEN];
  PetscTruth matoSet, matvSet;
  ierr = PetscOptionsGetString(PETSC_NULL, "-mato", matf, PETSC_MAX_PATH_LEN, &matoSet); 
           CHKERRQ(ierr);
  if (matoSet == PETSC_FALSE) {// put default name in matf; perhaps user set "-matv" only
    strcpy(matf, "pism_views");
  }
  strcpy(matlabOutVars, "\0");
  ierr = PetscOptionsGetString(PETSC_NULL, "-matv", matlabOutVars, PETSC_MAX_PATH_LEN, &matvSet); 
            CHKERRQ(ierr);
  if (matvSet == PETSC_TRUE) {
    strcat(matf, ".m");
    ierr = verbPrintf(1, grid.com, 
       "\n ... writing variables %s to Matlab file `%s'", matlabOutVars, matf); CHKERRQ(ierr);
    ierr = writeMatlabVars(matf); CHKERRQ(ierr); // see iMmatlab.cc
  }

  return 0;
}


PetscErrorCode IceModel::dumpToFile_netCDF(const char *fname) {
  PetscErrorCode ierr;

// bring in the result of applying pism-state-ncgen.py to pism_state.cdl
//   (see directory pism/src/netcdf/)
#include "../netcdf/write_attributes.c"

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.Mbz}; // Count local block: bed

  // Allocate some memory.
  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * PetscMax(grid.Mz, grid.Mbz);
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(double), &a_mpi); CHKERRQ(ierr);

// complete the output file
#include "../netcdf/complete_dump.cc"

  // We are done with these buffers
  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  return 0;
}


PetscErrorCode IceModel::dumpToFile_diagnostic_netCDF(const char *diag_fname) {
  PetscErrorCode ierr;

// bring in the result of applying pism-state-ncgen.py to pism_state.cdl
//   and pism_diag.fragment
//   (see directory pism/src/netcdf/)
// ("ncid" defined in included file)
#include "../netcdf/write_diag_attributes.c"

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.Mbz}; // Count local block: bed

  // Allocate some memory.
  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * PetscMax(grid.Mz, grid.Mbz);
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(double), &a_mpi); CHKERRQ(ierr);

// complete the output file as in usual dump
// ("uvel_id", "vvel_id", "wvel_id" defined in included file)
#include "../netcdf/complete_dump.cc"

  // now write additional 3-D diagnostic quantities
  ierr = u3.setVaridNC(uvel_id); CHKERRQ(ierr);
  ierr = u3.putVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = v3.setVaridNC(vvel_id); CHKERRQ(ierr);
  ierr = v3.putVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = w3.setVaridNC(wvel_id); CHKERRQ(ierr);
  ierr = w3.putVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  
  // We are done with these buffers
  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  return 0;
}


//! When reading a saved PISM model state, warn the user if options <tt>-Mx,-My,-Mz,-Mbz</tt> have been ignored.
PetscErrorCode IceModel::warnUserOptionsIgnored(const char *fname) {
  PetscErrorCode ierr;
  PetscInt       ignor;
  PetscTruth     M_Set;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mx", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -Mx ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-My", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -My ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mz", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -Mz ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mbz", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
              "WARNING: user option -Mbz ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  return 0;
}



//! Read a saved PISM model state in NetCDF format, for complete initialization of an evolution or diagnostic run.
/*! 
When initializing from a NetCDF input file, the input file determines 
the number of grid points (Mx,My,Mz,Mbz) and the dimensions (Lx,Ly,Lz) of the computational box.   
The user is warned when their command line options "-Mx", "-My", "-Mz", "-Mbz" are overridden.  
 */
PetscErrorCode IceModel::initFromFile_netCDF(const char *fname) {
  PetscErrorCode  ierr;
  int         ncid, stat;

  if (hasSuffix(fname, ".pb") == true) {
    SETERRQ1(1,"ERROR: .pb format not supported; cannot initialize from file %s", fname);
  }
  
  if (hasSuffix(fname, ".nc") == false) {
    ierr = verbPrintf(1,grid.com,
       "WARNING:  Unknown file format for %s.  Trying to read as NetCDF.\n",fname); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2,grid.com,"initializing from NetCDF file  %s  ...\n",
                     fname); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  size_t      dim[5];
  double      bdy[7];
  // note user option setting of -Lx,-Ly,-Lz will overwrite the corresponding settings from 
  //   this file but that the file's settings of Mx,My,Mz,Mbz will overwrite the user options 
  ierr = nct.get_dims_limits_lengths(ncid, dim, bdy, grid.com); CHKERRQ(ierr);
  grid.year = bdy[0] / secpera;
  grid.Mx = dim[1];
  grid.My = dim[2];
  grid.Mz = dim[3];
  grid.Mbz = dim[4];
  // grid.Lx, grid.Ly set from bdy[1], bdy[3] below in call to grid.rescale_using_zlevels()

  double *zlevs, *zblevs;
  zlevs = new double[grid.Mz];
  zblevs = new double[grid.Mbz];
  ierr = nct.get_vertical_dims(ncid, grid.Mz, grid.Mbz, zlevs, zblevs, grid.com); CHKERRQ(ierr);

  // re-allocate and fill grid.zlevels & zblevels with read values
  delete [] grid.zlevels;
  delete [] grid.zblevels;
  grid.zlevels = new PetscScalar[grid.Mz];
  grid.zblevels = new PetscScalar[grid.Mbz];
  for (PetscInt k = 0; k < grid.Mz; k++) {
    grid.zlevels[k] = (PetscScalar) zlevs[k];
  }
  for (PetscInt k = 0; k < grid.Mbz; k++) {
    grid.zblevels[k] = (PetscScalar) zblevs[k];
  }
  delete [] zlevs;  delete [] zblevs;  // done with these
  
  // set DA and create vecs; we have enough already to do so
  ierr = warnUserOptionsIgnored(fname); CHKERRQ(ierr);  
  ierr = grid.createDA(); CHKERRQ(ierr);
  // FIXME: note we *can* determine from the input file whether the hor. dims are truely periodic,
  // but this has not been done; here we simply require it is not periodic
  ierr = grid.rescale_using_zlevels(-bdy[1], -bdy[3]); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  
  // set IceModel::startYear, IceModel::endYear, grid.year, but respecting grid.year
  // which came from -if file, _unless_ -ys set by user
  ierr = setStartRunEndYearsFromOptions(PETSC_TRUE);  CHKERRQ(ierr);

  // Time to compute what we need.
  int s[] = {dim[0] - 1, grid.xs, grid.ys, 0};   // Start local block: t dependent; 
                                                 //   dim[0] is the number of t vals in file
  int c[] = {1, grid.xm, grid.ym, grid.Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.Mbz}; // Count local block: bed

  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.Mz;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(double), &a_mpi); CHKERRQ(ierr);

  // 2-D mapping
  ierr = nct.get_local_var(&grid, ncid, "lon", grid.da2, vLongitude, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "lat", grid.da2, vLatitude, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // 2-D model quantities: discrete
  ierr = nct.get_local_var(&grid, ncid, "mask", grid.da2, vMask, g2, 
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // 2-D model quantities: double
  ierr = nct.get_local_var(&grid, ncid, "usurf", grid.da2, vh, g2,  // DEPRECATED: pism_intent = diagnostic;
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);             //   WE SHOULD NOT BE READING THIS
                                                                              //   (CHOICES ABOUT WHAT -bif DOES ARE A
                                                                              //   DIFFERENT ISSUE)
  ierr = nct.get_local_var(&grid, ncid, "thk", grid.da2, vH, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "bwat", grid.da2, vHmelt, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "topg", grid.da2, vbed, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "dbdt", grid.da2, vuplift, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // Read vubarSSA and vvbarSSA if SSA is on, if not asked to ignore them and
  // if they are present in the input file.
  PetscTruth dontreadSSAvels = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL, "-dontreadSSAvels", &dontreadSSAvels); CHKERRQ(ierr);

  if ((grid.rank == 0) && useSSAVelocity && (!dontreadSSAvels)) {

    int flag;
    stat = nc_get_att_int(ncid, NC_GLOBAL, "haveSSAvelocities", &flag);
    if (stat == NC_NOERR)
      haveSSAvelocities = flag;
  }
  ierr = MPI_Bcast(&haveSSAvelocities, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  
  if (haveSSAvelocities == 1) {
    ierr = verbPrintf(2,grid.com,"Reading vubarSSA and vvbarSSA...\n"); CHKERRQ(ierr);

    ierr = nct.get_local_var(&grid, ncid, "vubarSSA", grid.da2, vubarSSA, g2,
			     s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
    ierr = nct.get_local_var(&grid, ncid, "vvbarSSA", grid.da2, vvbarSSA, g2,
			     s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  }

  // 2-D climate/bdry quantities
  ierr = nct.get_local_var(&grid, ncid, "artm", grid.da2, vTs, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "bheatflx", grid.da2, vGhf, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "acab", grid.da2, vAccum, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "tillphi", grid.da2, vtillphi, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);

  // 3-D model quantities
  ierr = T3.getVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = Tb3.getVecNC(ncid, s, cb, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = tau3.getVecNC(ncid, s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);

  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  // actually read the polar_stereographic if present
  if (grid.rank == 0) {
    int psid;
    stat = nc_inq_varid(ncid, "polar_stereographic", &psid); 
    if (stat == NC_NOERR) { // polar_stereo exists
      stat = nc_get_att_double(ncid, psid, "straight_vertical_longitude_from_pole",
                              &psParams.svlfp); CHKERRQ(nc_check(stat));
      stat = nc_get_att_double(ncid, psid, "latitude_of_projection_origin",
                              &psParams.lopo); CHKERRQ(nc_check(stat));
      stat = nc_get_att_double(ncid, psid, "standard_parallel",
                              &psParams.sp); CHKERRQ(nc_check(stat));
    }
  }
  // now broadcast whatever vals are on proc 0; this might be the default, which is fine
  ierr = MPI_Bcast(&psParams.svlfp, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&psParams.lopo, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&psParams.sp, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);

  // Get the current history length
  unsigned int history_len;	// used for communication
  if (grid.rank == 0) {
    size_t H;
    stat = nc_inq_attlen(ncid, NC_GLOBAL, "history", &H);
    history_len = (int)H;
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  // Broadcast the history length
  MPI_Bcast(&history_len, 1, MPI_INT, 0, grid.com);

  // Allocate some memory (if necessary)
  if (history_len > history_size - 1) {
    history_size = history_len + 1;
    delete[] history;
    history = new char[history_size];
  }

  // Zero out the allocated memory so that we don't need to worry about the
  // trailing zero in the history string (which might be absent).
  memset(history, 0, history_size);

  // Read the history string and close the file
  if (grid.rank == 0) {
    stat = nc_get_att_text(ncid, NC_GLOBAL, "history", history);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  // Broadcast the string
  MPI_Bcast(history, history_size, MPI_CHAR, 0, grid.com);

  initialized_p = PETSC_TRUE;
  return 0;
}


//! Manage regridding based on user options.  Call NCTool::regrid_local_var() to do each selected variable.
/*!
For each variable selected by option <tt>-regrid_vars</tt>, we regrid it onto the current grid from 
the NetCDF file specified by <tt>-regrid</tt>.

The default, if <tt>-regrid_vars</tt> is not given, is to regrid the 3 dimensional 
quantities \c tau3, \c T3, \c Tb3.  This is consistent with one standard purpose of 
regridding, which is to stick with current geometry through the downscaling procedure.  
Most of the time the user should carefully specify which variables to regrid.
 */
PetscErrorCode IceModel::regrid_netCDF(const char *regridFile) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];
  
  const int  npossible = 7;
  const char possible[20] = "bBehHLT";
  
  if (!hasSuffix(regridFile, ".nc")) {
    SETERRQ(1,"regridding is only possible if the source file has NetCDF format and extension '.nc'");
  }

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    strcpy(regridVars, "");
  }
  ierr = verbPrintf(2,grid.com, 
           "regridding variables with single character flags `%s' from NetCDF file `%s':", 
           regridVars,regridFile); CHKERRQ(ierr);

  // following are dimensions, limits and lengths, and id for *source* NetCDF file (regridFile)
  size_t dim[5];
  double bdy[7];
  int ncid, stat;

  // create "local interpolation context" from dimensions, limits, and lengths extracted from regridFile,
  //   and from information about the part of the grid owned by this processor
  if (grid.rank == 0) {
    stat = nc_open(regridFile, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  ierr = nct.get_dims_limits_lengths(ncid, dim, bdy, grid.com); CHKERRQ(ierr);  // see nc_util.cc
  // from regridFile: Mz = dim[3], Mbz = dim[4]  
  double *zlevs, *zblevs;
  zlevs = new double[dim[3]];
  zblevs = new double[dim[4]];
  ierr = nct.get_vertical_dims(ncid, dim[3], dim[4], zlevs, zblevs, grid.com); CHKERRQ(ierr);

  { // explicit scoping means destructor will be called for lic
    LocalInterpCtx lic(ncid, dim, bdy, zlevs, zblevs, grid);
    // ierr = lic.printGrid(grid.com); CHKERRQ(ierr);
    
    for (PetscInt k = 0; k < npossible; k++) {
      if (strchr(regridVars, possible[k])) {
       switch (possible[k]) {
         case 'b':
           ierr = verbPrintf(2, grid.com, "\n   b: regridding 'topg' ... "); CHKERRQ(ierr);
           ierr = nct.regrid_local_var("topg", 2, lic, grid, grid.da2, vbed, g2, false); CHKERRQ(ierr);
           break;
         case 'B':
           ierr = verbPrintf(2, grid.com, "\n   B: regridding 'litho_temp' ... "); CHKERRQ(ierr);
           ierr = Tb3.regridVecNC(4, lic);  CHKERRQ(ierr);
           break;
         case 'e':
           ierr = verbPrintf(2, grid.com, "\n   e: regridding 'age' ... "); CHKERRQ(ierr);
           ierr = tau3.regridVecNC(3, lic);  CHKERRQ(ierr);
           break;
         case 'h':
           ierr = verbPrintf(2, grid.com, "\n   h: regridding 'usurf' ... "); CHKERRQ(ierr);
           ierr = nct.regrid_local_var("usurf", 2, lic, grid, grid.da2, vh, g2, false); CHKERRQ(ierr);
           break;
         case 'H':
           ierr = verbPrintf(2, grid.com, "\n   H: regridding 'thk' ... "); CHKERRQ(ierr);
           ierr = nct.regrid_local_var("thk", 2, lic, grid, grid.da2, vH, g2, false); CHKERRQ(ierr);
           break;
         case 'L':
           ierr = verbPrintf(2, grid.com, "\n   L: regridding 'bwat' ... "); CHKERRQ(ierr);
           ierr = nct.regrid_local_var("bwat", 2, lic, grid, grid.da2, vHmelt, g2, false); CHKERRQ(ierr);
           break;
         case 'T':
           ierr = verbPrintf(2, grid.com, "\n   T: regridding 'temp' ... "); CHKERRQ(ierr);
           ierr = T3.regridVecNC(3, lic);  CHKERRQ(ierr);
           break;
       }
      }
    }

  }
  
  delete [] zlevs;  delete [] zblevs;

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  ierr = verbPrintf(2,grid.com, "\n"); CHKERRQ(ierr);
  return 0;
}

