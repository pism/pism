// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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

#include <petscsys.h>
#include <petscfix.h>


PetscErrorCode IceModel::regrid(const char *regridFile) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];
  InterpCtx i2, i3, i3b;

  ierr = PetscPrintf(grid.com, "regridding data from %s\n", regridFile); CHKERRQ(ierr);

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
  ierr = vPetscPrintf(grid.com, "  regrid interpolation context i2 created ... "); CHKERRQ(ierr);
  ierr = getInterpCtx(m.grid.da3, grid.da3, m, i3); CHKERRQ(ierr);
  ierr = vPetscPrintf(grid.com, "i3 created ... "); CHKERRQ(ierr);
  ierr = getInterpCtx(m.grid.da3b, grid.da3b, m, i3b); CHKERRQ(ierr);
  ierr = vPetscPrintf(grid.com, "i3b created.\n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    // As a default, we only regrid the mask and 3 dimensional quantities.  This
    // is consistent with the standard purpose which is to stick with current
    // geometry through the downscaling procedure.
    strcpy(regridVars, "mTBe");
  }

  ierr = regridVar(regridVars, 'm', i2, m.vMask, vMask);
  ierr = regridVar(regridVars, 'h', i2, m.vh, vh); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'H', i2, m.vH, vH); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'b', i2, m.vbed, vbed); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'a', i2, m.vAccum, vAccum); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 's', i2, m.vTs, vTs); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'g', i2, m.vGhf, vGhf); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'T', i3, m.vT, vT); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'B', i3b, m.vTb, vTb); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'e', i3, m.vtau, vtau); CHKERRQ(ierr);

  ierr = vPetscPrintf(grid.com, "\n"); CHKERRQ(ierr);

  ierr = destroyInterpCtx(i2); CHKERRQ(ierr);
  ierr = destroyInterpCtx(i3); CHKERRQ(ierr);
  ierr = destroyInterpCtx(i3b); CHKERRQ(ierr);

  PetscScalar run_years = endYear - startYear;
  ierr = setStartYear(m.grid.p->year); CHKERRQ(ierr);
  ierr = setRunYears(run_years); CHKERRQ(ierr);
  grid.p->year = startYear;
  
  // Stamp history
  ierr = stampHistoryString("<---- ---- Begin regrid block ---- ----"); CHKERRQ(ierr);
  ierr = stampHistoryString(m.grid.p->history); CHKERRQ(ierr);
  ierr = stampHistoryString(" ---- ----  End regrid block  ---- ---->\n"); CHKERRQ(ierr);
  
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
    } else if (Mzf == grid.p->Mbz || Mzf == 1) {
      ic.gc = cmodel.g3b; ic.gf = g3b;
      // If there is no basal component (Mzf == 1) then this actually doesn't
      // matter since it gets overwritten by the boundary condition in the
      // temperature equation, so we will just do something.
      zl = -grid.p->Lbz; zh = 0.0;
      zzl = -cgrid.p->Lbz; zzh = 0.0;
      if (zl == 0.0) zl = -1.0;
      if (zzl == 0.0) zzl = -1.0;
    } else {
      SETERRQ(1, "Something is broken here.");
    }
  }

  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE,
                         Mxf * Myf * Mzf, Mxc * Myc * Mzc,
                         8, PETSC_NULL, 8, PETSC_NULL, &(ic.A)); CHKERRQ(ierr);
  ierr = MatGetVecs(ic.A, &(ic.coarse), &(ic.fine)); CHKERRQ(ierr);

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
  
  // This is a bit ugly.  We need to move the data on local vectors to a
  // global vector which is compatible with the interpolation matrix `ic.A'.
  // Then we move this onto a local vector.  We should be able to do this
  // cleanly with scatters, but my first effort was repelled.  If performance
  // is an issue, then we should improve this part.  Since regridding is only
  // done once per model run, it should not be a problem considering that we
  // do it in parallel.
  ierr = vPetscPrintf(grid.com, "  regridding %c ... ",c); CHKERRQ(ierr);

  ierr = DALocalToGlobal(ic.dac, src, INSERT_VALUES, ic.gc); CHKERRQ(ierr);

  PetscInt low, high, n;
  PetscInt *ida;
  PetscScalar *a;
  PetscInt xs, ys, zs, xm, ym, zm, Mz, My, Mx;
  ierr = DAGetCorners(ic.dac, &zs, &ys, &xs, &zm, &ym, &xm); CHKERRQ(ierr);
  ierr = DAGetInfo(ic.dac, PETSC_NULL, &Mz, &My, &Mx, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  ida = new PetscInt[xm * ym * zm];
  n = 0;
  for (PetscInt i = xs; i < xs + xm; i++) {
    for (PetscInt j = ys; j < ys + ym; j++) {
      for (PetscInt k = zs; k < zs + zm; k++) {
        ida[n++] = i * (My * Mz) + j * Mz + k;
      }
    }
  }
    
  ierr = VecGetArray(ic.gc, &a); CHKERRQ(ierr);
  ierr = VecSetValues(ic.coarse, n, ida, a, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecRestoreArray(ic.gc, &a); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ic.coarse); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ic.coarse); CHKERRQ(ierr);
  delete [] ida;
    
  // ierr = VecView(gsrc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // ierr = PetscPrintf(grid.com, "#######################################\n"); CHKERRQ(ierr);
  // ierr = VecView(interpCtx.coarse, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = MatMult(ic.A, ic.coarse, ic.fine); CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(ic.fine, &low, &high); CHKERRQ(ierr);
  ierr = VecGetArray(ic.fine, &a); CHKERRQ(ierr);
  ida = new PetscInt[high - low];
  for (PetscInt i = 0; i < high - low; i++) {
    ida[i] = low + i;
  }
  ierr = VecSetValues(ic.gf, high - low, ida, a, INSERT_VALUES);
  ierr = VecAssemblyBegin(ic.gf); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ic.gf); CHKERRQ(ierr);
  delete [] ida;
  ierr = VecRestoreArray(ic.fine, &a); CHKERRQ(ierr);

  ierr = DAGlobalToLocalBegin(ic.daf, ic.gf, INSERT_VALUES, dest); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(ic.daf, ic.gf, INSERT_VALUES, dest); CHKERRQ(ierr);

  return 0;
}

