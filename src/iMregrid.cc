// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cstring>
#include "iceModel.hh"
#include "nc_util.hh"

#include <petscsys.h>
#include <petscfix.h>

// The regrid process has been cleaned up, but the fundamentals are the same as
// before.  We need to move a local vector from a coarse grid to a fine grid
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

// Possible bugs:
//   VecGetLocalSize() is not guaranteed to work on all platforms
//   MatGetOwnershipRange() is not guaranteed to work on all platforms


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

  //ierr = get_dimensions(ncid, dim, bdy, grid.com); CHKERRQ(ierr);
  ierr = get_dimensions(ncid, dim, bdy, &bdy_time, grid.com); CHKERRQ(ierr);

  // Get Local Interpolation Context
  LocalInterpCtx lic;
  ierr = get_LocalInterpCtx(ncid, dim, bdy, bdy_time, lic, grid); CHKERRQ(ierr);

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

  ierr = verbPrintf(2,grid.com, "\n"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::regrid(const char *regridFile) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];
  InterpCtx i2, i3, i3b;

  if (hasSuffix(regridFile, ".nc")) {
    // We have a much better method for regridding netCDF files, above
    ierr = regrid_netCDF(regridFile); CHKERRQ(ierr);
    return 0;
  }

  ierr = verbPrintf(2,grid.com, "regridding data from `%s'\n", regridFile); CHKERRQ(ierr);

  IceGrid g(grid.com, grid.rank, grid.size);
  IceModel m(g, ice);
  m.setShowViewers(PETSC_FALSE);
  m.setAllowRegridding(PETSC_FALSE);
  m.initFromFile(regridFile);
  m.afterInitHook();

  if ((grid.p->Lx > m.grid.p->Lx)
      || (grid.p->Ly > m.grid.p->Ly)
      || (grid.p->Lz > m.grid.p->Lz)
      || (grid.p->Lbz > m.grid.p->Lbz)) {
    SETERRQ(1, "New grid exceeds extent of old grid.");
  }

  ierr = getInterpCtx(m.grid.da2, grid.da2, m, i2); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, "  regrid interpolation context i2 created ... "); CHKERRQ(ierr);
  ierr = getInterpCtx(m.grid.da3, grid.da3, m, i3); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, "i3 created ... "); CHKERRQ(ierr);
  ierr = getInterpCtx(m.grid.da3b, grid.da3b, m, i3b); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, "i3b created.\n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    // As a default, we only regrid the 3 dimensional quantities.  This
    // is consistent with the standard purpose which is to stick with current
    // geometry through the downscaling procedure.
    strcpy(regridVars, "TBe");
  }

  ierr = regridVar(regridVars, 'm', i2, m.vMask, vMask);  // not recommended?
  // ierr = regridVar(regridVars, 'h', i2, m.vh, vh); CHKERRQ(ierr); // it is diagnositic 

  // regridable model state variables
  ierr = regridVar(regridVars, 'b', i2, m.vbed, vbed); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'B', i3b, m.vTb, vTb); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'e', i3, m.vtau, vtau); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'H', i2, m.vH, vH); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'L', i2, m.vHmelt, vHmelt); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'T', i3, m.vT, vT); CHKERRQ(ierr);

  // regriddable climate variables
  ierr = regridVar(regridVars, 'a', i2, m.vAccum, vAccum); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'g', i2, m.vGhf, vGhf); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 's', i2, m.vTs, vTs); CHKERRQ(ierr);

  ierr = verbPrintf(3, grid.com, "\n"); CHKERRQ(ierr);

  ierr = destroyInterpCtx(i2); CHKERRQ(ierr);
  ierr = destroyInterpCtx(i3); CHKERRQ(ierr);
  ierr = destroyInterpCtx(i3b); CHKERRQ(ierr);

  ierr = setStartYear(m.grid.p->year); CHKERRQ(ierr);
  ierr = setEndYear(endYear); CHKERRQ(ierr);
  grid.p->year = startYear;
  
  // Stamp history
  ierr = stampHistoryString("<---- ---- Begin regrid block ---- ----"); CHKERRQ(ierr);
  ierr = stampHistoryString(m.grid.p->history); CHKERRQ(ierr);
  ierr = stampHistoryString(" ---- ----  End regrid block  ---- ---->\n"); CHKERRQ(ierr);
  
  ierr = verbPrintf(2,grid.com, "regridding done\n"); CHKERRQ(ierr);

  // The model that we just grabbed data from will fall out of scope and its
  // memory will be freed at the end of this method.  
  return 0;
}


int ind(double i, double j, double k, int N, int K) {
  return (int)i * N * K + (int)j * K + (int)k;
}


PetscErrorCode IceModel::getInterpCtx(const DA dac, const DA daf,
                                      const IceModel &cmodel, InterpCtx &ic) {
  PetscErrorCode ierr;
  PetscInt Mzc, Myc, Mxc, Mzf, Myf, Mxf;
  PetscInt rowl, rowh;
  PetscReal xl, xh, xxl, xxh;
  PetscReal yl, yh, yyl, yyh;
  PetscReal zl, zh, zzl, zzh;
  
  ic.dac = dac;
  ic.daf = daf;
  
  ierr = DAGetInfo(dac, PETSC_NULL, &Mzc, &Myc, &Mxc, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DAGetInfo(daf, PETSC_NULL, &Mzf, &Myf, &Mxf, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  const IceGrid &cgrid = cmodel.grid;
  // The concept of x, y, z are goofy here so that we don't have to duplicate a
  // bunch of code.  We have to make sure that the span (*h - *l) is always positive
  if (Mxf == 1) { // We are working with a 2 dimensional array
    ic.gc = cmodel.g2; ic.gf = g2;
    zl = -grid.p->Ly; zh = grid.p->Ly;
    yl = -grid.p->Lx; yh = grid.p->Lx;
    xl = 0.0; xh = 1.0;
    zzl = -cgrid.p->Ly; zzh = cgrid.p->Ly;
    yyl = -cgrid.p->Lx; yyh = cgrid.p->Lx;
    xxl = 0.0; xxh = 1.0;
  } else { // It is 3 dimensional
    xl = -grid.p->Lx; xh = grid.p->Lx;
    yl = -grid.p->Ly; yh = grid.p->Ly;
    xxl = -cgrid.p->Lx; xxh = cgrid.p->Lx;
    yyl = -cgrid.p->Ly; yyh = cgrid.p->Ly;
    if (Mzf == grid.p->Mz) {
      ic.gc = cmodel.g3; ic.gf = g3;
      zl = 0.0; zh = grid.p->Lz;
      zzl = 0.0; zzh = cgrid.p->Lz;
    } else if (Mzf == grid.p->Mbz) {
      ic.gc = cmodel.g3b; ic.gf = g3b;
      zl = -grid.p->Lbz; zh = 0.0;
      zzl = -cgrid.p->Lbz; zzh = 0.0;
      if (zl == 0.0) zl = -1.0;
      if (zzl == 0.0) zzl = -1.0;
    } else {
      SETERRQ(1, "Something is broken here.");
    }
  }

  ierr = DACreateNaturalVector(ic.dac, &(ic.coarse)); CHKERRQ(ierr);
  ierr = DACreateNaturalVector(ic.daf, &(ic.fine)); CHKERRQ(ierr);
  PetscInt cLocSize, fLocSize;
  ierr = VecGetLocalSize(ic.coarse, &cLocSize); CHKERRQ(ierr);
  ierr = VecGetLocalSize(ic.fine, &fLocSize); CHKERRQ(ierr);
  ierr = MatCreateMPIAIJ(grid.com, fLocSize, cLocSize,
                         PETSC_DETERMINE, PETSC_DETERMINE,
                         8, PETSC_NULL, 8, PETSC_NULL, &(ic.A)); CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(ic.A, &rowl, &rowh); CHKERRQ(ierr);

  for (PetscInt row = rowl; row < rowh; row++) {
    const PetscInt i = row / (Myf * Mzf);
    const PetscInt j = (row - i * (Myf * Mzf)) / Mzf;
    const PetscInt k = row % Mzf;

    // We are multiplexing the 2D and 3D cases and here we do not want to divide
    // by zero in the 2D case.
    const PetscReal x = (Mxf == 1) ? 0 : xl + i * (xh - xl) / (Mxf - 1);
    const PetscReal y = yl + j * (yh - yl) / (Myf - 1);
    const PetscReal z = zl + k * (zh - zl) / (Mzf - 1);

    const PetscReal ii = (Mxc - 1) * (x - xxl) / (xxh - xxl);
    const PetscReal jj = (Myc - 1) * (y - yyl) / (yyh - yyl);
    const PetscReal kk = (Mzc - 1) * (z - zzl) / (zzh - zzl);

    // These are in [0,1)
    const PetscReal xx = ii - floor(ii);
    const PetscReal yy = jj - floor(jj);
    const PetscReal zz = kk - floor(kk);
    
    const PetscInt col[] = {
      ind(floor(ii), floor(jj), floor(kk), Myc, Mzc),
      ind(floor(ii), floor(jj), ceil(kk), Myc, Mzc),
      ind(floor(ii), ceil(jj), floor(kk), Myc, Mzc),
      ind(floor(ii), ceil(jj), ceil(kk), Myc, Mzc),
      ind(ceil(ii), floor(jj), floor(kk), Myc, Mzc),
      ind(ceil(ii), floor(jj), ceil(kk), Myc, Mzc),
      ind(ceil(ii), ceil(jj), floor(kk), Myc, Mzc),
      ind(ceil(ii), ceil(jj), ceil(kk), Myc, Mzc)
    };
    const PetscScalar weight[] = {
      (1 - xx) * (1 - yy) * (1 - zz),
      (1 - xx) * (1 - yy) * zz,
      (1 - xx) * yy * (1 - zz),
      (1 - xx) * yy * zz,
      xx * (1 - yy) * (1 - zz),
      xx * (1 - yy) * zz,
      xx * yy * (1 - zz),
      xx * yy * zz
    };

    ierr = MatSetValues(ic.A, 1, &row, 8, col, weight, ADD_VALUES); CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(ic.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ic.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::destroyInterpCtx(InterpCtx &i) {
  PetscErrorCode ierr;

  ierr = MatDestroy(i.A); CHKERRQ(ierr);
  ierr = VecDestroy(i.coarse); CHKERRQ(ierr);
  ierr = VecDestroy(i.fine); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::regridVar(const char *vars, char c, const InterpCtx &ic,
                                   const Vec src, Vec dest) {
  PetscErrorCode ierr;

  if (! strchr(vars, c)) {
    return 0;
  }
  
  ierr = verbPrintf(3, grid.com, "  regridding %c ...\n",c); CHKERRQ(ierr);

  ierr = DALocalToGlobal(ic.dac, src, INSERT_VALUES, ic.gc); CHKERRQ(ierr);

  // ierr = VecView(gsrc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // ierr = PetscPrintf(grid.com, "#######################################\n"); CHKERRQ(ierr);
  // ierr = VecView(interpCtx.coarse, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = DAGlobalToNaturalBegin(ic.dac, ic.gc, INSERT_VALUES, ic.coarse); CHKERRQ(ierr);
  ierr = DAGlobalToNaturalEnd(ic.dac, ic.gc, INSERT_VALUES, ic.coarse); CHKERRQ(ierr);

  PetscInt mm, mn, fm, cn;
  ierr = MatGetSize(ic.A, &mm, &mn); CHKERRQ(ierr);
  ierr = VecGetSize(ic.fine, &fm); CHKERRQ(ierr);
  ierr = VecGetSize(ic.coarse, &cn); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, "Interpolation matrix size: %d, %d\n", mm, mn); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, "Coarse vector size: %d\n", cn); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, "Fine vector size: %d\n", fm); CHKERRQ(ierr);
  
  ierr = MatMult(ic.A, ic.coarse, ic.fine); CHKERRQ(ierr);

  ierr = DANaturalToGlobalBegin(ic.daf, ic.fine, INSERT_VALUES, ic.gf); CHKERRQ(ierr);
  ierr = DANaturalToGlobalEnd(ic.daf, ic.fine, INSERT_VALUES, ic.gf); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(ic.daf, ic.gf, INSERT_VALUES, dest); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(ic.daf, ic.gf, INSERT_VALUES, dest); CHKERRQ(ierr);

  return 0;
}

