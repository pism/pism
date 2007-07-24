// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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
#include "nc_util.hh"
#include "iceModel.hh"

// We need to move a local vector from a coarse grid to a fine grid
// (there is no requirement that the `coarse' grid need actually be coarser than
// the `fine' grid).  Things are really ugly in any ordering other than the
// `natural' ordering, so we move local -> global -> natural with the coarse
// data.  We must have defined a weighting matrix which operates on this coarse
// natural vector to produce a fine natural vector.  Currently, we use
// tri-linear interpolation.  After applying the matrix, we move the fine vector
// back to a local vector: natural -> global -> local.  It is theoretically
// possible to make the matrix operate on vectors in the Petsc global ordering,
// but that seems like a mess.  In particular, since the matrix is not square,
// we cannot use DAGetMatrix() or the like.

PetscErrorCode IceModel::regrid(const char *regridFile) {
  PetscErrorCode ierr;

  if (hasSuffix(regridFile, ".nc")) {
    ierr = regrid_netCDF(regridFile); CHKERRQ(ierr);
    return 0;
  } else {
    SETERRQ(1,"unable to regrid non-NetCDF file");
  }
}
  
PetscErrorCode IceModel::regrid_netCDF(const char *regridFile) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(2,grid.com, "regridding data from `%s': ", regridFile); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    // As a default, we only regrid the 3 dimensional quantities.  This
    // is consistent with one standard purpose which is to stick with current
    // geometry through the downscaling procedure.
    strcpy(regridVars, "TBe");
  }

  size_t dim[5];
  float bdy[7];
  double bdy_time;
  int ncid, stat;

  if (grid.rank == 0) {
    stat = nc_open(regridFile, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  ierr = get_dimensions(ncid, dim, bdy, &bdy_time, grid.com); CHKERRQ(ierr);  // see nc_util.cc

  // Get Local Interpolation Context
  LocalInterpCtx lic;
  ierr = get_LocalInterpCtx(ncid, dim, bdy, bdy_time, lic, grid); CHKERRQ(ierr);  // see nc_util.cc

  ierr = regrid_local_var(regridVars, 'h', "h", 2, lic, grid, grid.da2, vh, g2);  // see nc_util.cc
  CHKERRQ(ierr);
  ierr = regrid_local_var(regridVars, 'H', "H", 2, lic, grid, grid.da2, vH, g2);
  CHKERRQ(ierr);
  ierr = regrid_local_var(regridVars, 'L', "Hmelt", 2, lic, grid, grid.da2, vHmelt, g2);
  CHKERRQ(ierr);
  ierr = regrid_local_var(regridVars, 'b', "b", 2, lic, grid, grid.da2, vbed, g2);
  CHKERRQ(ierr);
  ierr = regrid_local_var(regridVars, 'T', "T", 3, lic, grid, grid.da3, vT, g3);
  CHKERRQ(ierr);
  ierr = regrid_local_var(regridVars, 'B', "Tb", 4, lic, grid, grid.da3b, vTb, g3b);
  CHKERRQ(ierr);
  ierr = regrid_local_var(regridVars, 'e', "age", 3, lic, grid, grid.da3, vtau, g3);
  CHKERRQ(ierr);

  ierr = PetscFree(lic.a); CHKERRQ(ierr);
  
  if (grid.rank == 0) {
    // If we want history from the regridded file, we can get it here and broadcast it.
    // stat = nc_get_att_text(ncid, NC_GLOBAL, "history", grid.p->history);
    // CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  ierr = verbPrintf(2,grid.com, "\n"); CHKERRQ(ierr);

  return 0;
}

