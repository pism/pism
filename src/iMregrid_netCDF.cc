// Copyright (C) 2007 Jed Brown
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#if (WITH_NETCDF)

#include <cstring>
#include "iceModel.hh"
#include "nc_util.hh"

PetscErrorCode IceModel::regrid_netCDF(const char *regridFile) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    // As a default, we only regrid the 3 dimensional quantities.  This
    // is consistent with the standard purpose which is to stick with current
    // geometry through the downscaling procedure.
    strcpy(regridVars, "TBe");
  }

  size_t dim[5];
  float bdy[7];
  int ncid, stat;

  if (grid.rank == 0) {
    stat = nc_open(regridFile, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  ierr = get_dimensions(ncid, dim, bdy, grid.com); CHKERRQ(ierr);

  // Get Local Interpolation Context
  LocalInterpCtx lic;
  ierr = get_LocalInterpCtx(ncid, dim, bdy, lic, grid); CHKERRQ(ierr);

  ierr = regrid_local_var(regridVars, 'h', "h", 2, lic, grid, grid.da2, vh, g2);
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

  return 0;
}

#endif
