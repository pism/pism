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


//! Build the VecScatter and indexing scheme, from processor zero to others, used in NetCDF IO.
/*!
This procedure is used for NetCDF input and output but does not actually refer to any 
NetCDF constructs.  Only addresses the scatter from processor zero.  It could, 
potentially, be used by other code which 
does the same scatter to processor zero.
 */
PetscErrorCode IceModel::getIndZero(DA da, Vec vind, Vec vindzero, VecScatter ctx) {
  PetscErrorCode ierr;
  PetscInt low, high, n;
  PetscInt *ida;
  PetscScalar *a;
  PetscInt xs, ys, xm, ym, My, Mx;
  ierr = DAGetCorners(da, &ys, &xs, PETSC_NULL, &ym, &xm, PETSC_NULL); CHKERRQ(ierr);
  ierr = DAGetInfo(da, PETSC_NULL, &My, &Mx, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(vind, &low, &high); CHKERRQ(ierr);

  ida = new PetscInt[xm * ym];
  a = new PetscScalar[xm * ym];
  n = 0;
  for (PetscInt i = xs; i < xs + xm; i++) {
    for (PetscInt j = ys; j < ys + ym; j++) {
      ida[n] = i * My + j;
      a[n] = low + (i - xs) * ym + (j - ys);
      n++;
    }
  }
  if (n != high - low) {
    SETERRQ(1, "This should not happen");
  }
    
  ierr = VecSetValues(vind, n, ida, a, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vind); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vind); CHKERRQ(ierr);

  delete [] ida;
  delete [] a;
  ierr = VecScatterBegin(ctx, vind, vindzero, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx, vind, vindzero, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  return 0;
}


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
  // note grid.p->year has already been set from input file
  PetscScalar usrStartYear, usrEndYear, usrRunYears;
  PetscTruth ysSet, yeSet, ySet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ys", &usrStartYear, &ysSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ye", &usrEndYear, &yeSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-y", &usrRunYears, &ySet); CHKERRQ(ierr);
  if (ysSet == PETSC_TRUE) {
    // user option overwrites data
    ierr = setStartYear(usrStartYear); CHKERRQ(ierr);
    grid.p->year = usrStartYear;
  } else if (grid_p_year_VALID == PETSC_TRUE) {
    ierr = setStartYear(grid.p->year); CHKERRQ(ierr);
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
    ierr = setEndYear(DEFAULT_RUN_YEARS + startYear); CHKERRQ(ierr);
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
      ierr = verbPrintf(1, grid.com, 
            "Writing model state, with full 3D velocities, to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_diagnostic_netCDF(ncf); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(1, grid.com, "Writing model state to file `%s'", ncf); CHKERRQ(ierr);
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

// bring in the result of applying ncgen.rb to pism_state.cdl (see directory pism/src/netcdf/)
#include "../netcdf/write_attributes.c"

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.p->Mbz}; // Count local block: bed

//  // Allocate some memory.  We will assume that 3d ice vectors (Mx x My x Mz) are the largest.
  // Allocate some memory.
  void *a_mpi;
  int a_len, max_a_len;
//  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
  max_a_len = a_len = grid.xm * grid.ym * PetscMax(grid.p->Mz, grid.p->Mbz);
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

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

// bring in the result of applying ncgen.rb to pism_state.cdl along with pism_diag.fragment
// ("ncid" defined in included file)
#include "../netcdf/write_diag_attributes.c"

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.p->Mbz}; // Count local block: bed

//  // Allocate some memory.  We will assume that vectors based on grid.da3 are the largest.
  // Allocate some memory.
  void *a_mpi;
  int a_len, max_a_len;
//  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
  max_a_len = a_len = grid.xm * grid.ym * PetscMax(grid.p->Mz, grid.p->Mbz);
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

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
  size_t      dim[5];
  float       bdy[7];
  double 	  bdy_time;
  int         ncid, stat;

  if (hasSuffix(fname, ".pb") == true) {
    SETERRQ1(1,"ERROR: .pb format not supported; cannot initialize from file %s", fname);
  }
  
  if (hasSuffix(fname, ".nc") == false) {
    ierr = verbPrintf(1,grid.com,
       "WARNING:  Unknown file format for %s.  Trying to read as NetCDF.\n",fname); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2,grid.com,"initializing from NetCDF format file  %s  ...\n",
                     fname); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  // note user option setting of -Lx,-Ly,-Lz will overwrite the corresponding settings from 
  //   this file but that the file's settings of Mx,My,Mz,Mbz will overwrite the user options 
  ierr = nct.get_dims_limits_lengths(ncid, dim, bdy, &bdy_time, grid.com); CHKERRQ(ierr);
  // grid.p->year = bdy[0] / secpera;  // this was the float version; bdy_time is the double version
  grid.p->year = bdy_time / secpera;
  grid.p->Mx = dim[1];
  grid.p->My = dim[2];
  grid.p->Mz = dim[3];
  grid.p->Mbz = dim[4];
  
  float *zlevs, *zblevs;
  zlevs = new float[grid.p->Mz];
  zblevs = new float[grid.p->Mbz];
  ierr = nct.get_vertical_dims(ncid, grid.p->Mz, grid.p->Mbz, zlevs, zblevs, grid.com); CHKERRQ(ierr);

  // re-allocate and fill grid.zlevels & zblevels with read values
  delete [] grid.zlevels;
  delete [] grid.zblevels;
  grid.zlevels = new PetscScalar[grid.p->Mz];
  grid.zblevels = new PetscScalar[grid.p->Mbz];
  for (PetscInt k = 0; k < grid.p->Mz; k++) {
    grid.zlevels[k] = (PetscScalar) zlevs[k];
  }
  for (PetscInt k = 0; k < grid.p->Mbz; k++) {
    grid.zblevels[k] = (PetscScalar) zblevs[k];
  }
  delete [] zlevs;  delete [] zblevs;  // done with these
  
  // set DA and create vecs; we have enough already to do so
  ierr = warnUserOptionsIgnored(fname); CHKERRQ(ierr);  
  ierr = grid.createDA(); CHKERRQ(ierr);
  // FIXME: note we *can* determine from the input file whether the hor. dims are truely periodic,
  // but this has not been done; here we simply require it is not periodic
  ierr = grid.rescale_using_zlevels(-bdy[1], -bdy[3], PETSC_FALSE); CHKERRQ(ierr);  // must go after createDA() ...
  ierr = createVecs(); CHKERRQ(ierr);
  
  // set IceModel::startYear, IceModel::endYear, grid.p->year, but respecting grid.p->year
  // which came from -if file, _unless_ -ys set by user
  ierr = setStartRunEndYearsFromOptions(PETSC_TRUE);  CHKERRQ(ierr);

  // Time to compute what we need.
  int s[] = {dim[0] - 1, grid.xs, grid.ys, 0};   // Start local block: t dependent; dim[0] is the number of t vals in file
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.p->Mbz}; // Count local block: bed

  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

  // these are treated like 2-D constant quantities
  ierr = nct.get_local_var(&grid, ncid, "lon", NC_FLOAT, grid.da2, vLongitude, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "lat", NC_FLOAT, grid.da2, vLatitude, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D model quantities
  ierr = nct.get_local_var(&grid, ncid, "mask", NC_BYTE, grid.da2, vMask, g2, s, c, 3,
                       a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "usurf", NC_FLOAT, grid.da2, vh, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "thk", NC_FLOAT, grid.da2, vH, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "bwat", NC_FLOAT, grid.da2, vHmelt, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "topg", NC_FLOAT, grid.da2, vbed, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "dbdt", NC_FLOAT, grid.da2, vuplift, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D climate quantities
  ierr = nct.get_local_var(&grid, ncid, "artm", NC_FLOAT, grid.da2, vTs, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "bheatflx", NC_FLOAT, grid.da2, vGhf, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = nct.get_local_var(&grid, ncid, "acab", NC_FLOAT, grid.da2, vAccum, g2,
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

  // read the history *and* close the file
  int history_len;
  if (grid.rank == 0) {
    stat = nc_get_att_text(ncid, NC_GLOBAL, "history", grid.p->history);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
    history_len = strnlen(grid.p->history, HISTORY_STRING_LENGTH);
  }
  MPI_Bcast(&history_len, 1, MPI_INT, 0, grid.com);
  MPI_Bcast(grid.p->history, history_len, MPI_CHAR, 0, grid.com);

  setConstantGrainSize(DEFAULT_GRAIN_SIZE);

//  ierr = reconfigure_legacy_Mbz(); CHKERRQ(ierr);

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

  if (!hasSuffix(regridFile, ".nc")) {
    SETERRQ(1,"regridding is only possible if the source file has NetCDF format and extension '.nc'");
  }

  ierr = verbPrintf(1,grid.com,
          "WARNING: regrid_netCDF() CURRENTLY assumes the regrid file has evenly spaced vertical coordinate!\n");
          CHKERRQ(ierr); 

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    strcpy(regridVars, "TBe");
  }
  ierr = verbPrintf(2,grid.com, "regridding variables with single character flags `%s' from NetCDF file `%s'", 
                    regridVars,regridFile); CHKERRQ(ierr);

  // following are dimensions, limits and lengths, and id for *source* NetCDF file (regridFile)
  size_t dim[5];
  float bdy[7];
  double bdy_time;
  int ncid, stat;

  // create "local interpolation context" from dimensions, limits, and lengths extracted from regridFile,
  //   and from information about the part of the grid owned by this processor
  if (grid.rank == 0) {
    stat = nc_open(regridFile, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  ierr = nct.get_dims_limits_lengths(ncid, dim, bdy, &bdy_time, grid.com); CHKERRQ(ierr);  // see nc_util.cc
  // from regridFile: Mz = dim[3], Mbz = dim[4]  
  float *zlevs, *zblevs;
  zlevs = new float[dim[3]];
  zblevs = new float[dim[4]];
  ierr = nct.get_vertical_dims(ncid, dim[3], dim[4], zlevs, zblevs, grid.com); CHKERRQ(ierr);
  LocalInterpCtx lic;
  ierr = nct.form_LocalInterpCtx(ncid, dim, bdy, bdy_time, zlevs, zblevs, lic, grid); CHKERRQ(ierr);  // see nc_util.cc
  delete [] zlevs;  delete [] zblevs;

  // do each variable
  ierr = nct.regrid_local_var(regridVars, 'h', "usurf", 2, lic, grid, grid.da2, vh, g2);  // see nc_util.cc
  CHKERRQ(ierr);
  ierr = nct.regrid_local_var(regridVars, 'H', "thk", 2, lic, grid, grid.da2, vH, g2);
  CHKERRQ(ierr);
  ierr = nct.regrid_local_var(regridVars, 'L', "bwat", 2, lic, grid, grid.da2, vHmelt, g2);
  CHKERRQ(ierr);
  ierr = nct.regrid_local_var(regridVars, 'b', "topg", 2, lic, grid, grid.da2, vbed, g2);
  CHKERRQ(ierr);
  ierr = nct.regrid_local_var(regridVars, 'a', "acab", 2, lic, grid, grid.da2, vAccum, g2);
  CHKERRQ(ierr);
  ierr = T3.regridVecNC(  regridVars, 'T', 3, lic);  CHKERRQ(ierr);
  ierr = tau3.regridVecNC(regridVars, 'e', 3, lic);  CHKERRQ(ierr);
  ierr = Tb3.regridVecNC( regridVars, 'B', 4, lic);  CHKERRQ(ierr);

  ierr = PetscFree(lic.zlevs); CHKERRQ(ierr);  
  ierr = PetscFree(lic.zblevs); CHKERRQ(ierr);  
  ierr = PetscFree(lic.a); CHKERRQ(ierr);  
  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  ierr = verbPrintf(2,grid.com, "\n"); CHKERRQ(ierr);
  return 0;
}

