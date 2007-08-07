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
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"
#include "iceModel.hh"

PetscErrorCode IceModel::cleanJustNan_legacy() {
  // this procedure cleans data tied to Antarctica, including BAS data (acquired for PISM in 2004)
  // added in rev152 for ant_2007/ant_init.nc
  PetscErrorCode  ierr;
  PetscScalar     **h, **H, **bed, **ac, **Ts;

  /*
  * Clean nan values from h, H, bed, Ts, accum
  * Also initialize temperature at depth and set mask using floatation criterion.
  */
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &ac); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      /*
      * The undefined points should only be far out in the ocean. We probably don't
      * want ice there anyway so we provide values that will unobtrusively allow
      * transient ice shelves.
      */
      if (isnan(h[i][j])) h[i][j] = DEFAULT_h_VALUE_MISSING;
      if (isnan(H[i][j])) H[i][j] = DEFAULT_H_VALUE_MISSING;
      if (isnan(bed[i][j])) bed[i][j] = DEFAULT_BED_VALUE_MISSING;
      if (isnan(ac[i][j])) ac[i][j] = DEFAULT_ACCUM_VALUE_MISSING;
      if (isnan(Ts[i][j])) Ts[i][j] = DEFAULT_SURF_TEMP_VALUE_MISSING;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &ac); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::cleanInputData_legacy() {
  // this procedure cleans data tied to Antarctica, including BAS data (acquired for PISM in 2004)
  PetscErrorCode  ierr;
  PetscScalar     **h, **H, **bed, **ac, **Ts;

  // Convert temperature from Celsius to Kelvin
  ierr = VecShift(vTs, ice.meltingTemp); CHKERRQ(ierr);

  // Change units from m/a to m/s
  ierr = VecScale(vAccum, 1/secpera); CHKERRQ(ierr);
  ierr = VecScale(vuplift, 1/secpera); CHKERRQ(ierr);
  ierr = VecScale(vbalvel, 1/secpera); CHKERRQ(ierr);

  /*
  * Clean nan values from h, H, bed, Ts, accum
  * Also initialize temperature at depth and set mask using floatation criterion.
  */
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &ac); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      /*
      * The undefined points should only be far out in the ocean. We probably don't
      * want ice there anyway so we provide values that will unobtrusively allow
      * transient ice shelves.
      */
      if (isnan(h[i][j])) h[i][j] = DEFAULT_h_VALUE_MISSING;
      if (isnan(H[i][j])) H[i][j] = DEFAULT_H_VALUE_MISSING;
      if (isnan(bed[i][j])) bed[i][j] = DEFAULT_BED_VALUE_MISSING;
      if (isnan(ac[i][j])) ac[i][j] = DEFAULT_ACCUM_VALUE_MISSING;
      if (isnan(Ts[i][j])) Ts[i][j] = DEFAULT_SURF_TEMP_VALUE_MISSING;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &ac); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::createMask_legacy(PetscTruth balVelRule) {
  PetscErrorCode  ierr;
  Vec             vsliding;
  PetscScalar     **mask, **balvel, **sliding, **ubar, **vbar, **accum, **H;

  if (balVelRule == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com, 
                       "creating modal mask using balance velocities from input file ... "); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, 
                       "creating modal mask by simple floating/grounded decision ... "); CHKERRQ(ierr);
  }

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      /* The BAS input mask uses:
      *    0 : Ocean           1 : Grounded
      *    2 : Floating        3 : Rock outcrop
      * We lose information about "Grounded" vs "Rock outcrop"
      */
      switch (intMask(mask[i][j])) {
        case 0:
          mask[i][j] = MASK_FLOATING_OCEAN0;
          // FIXME: this ablation mechanism should be replaced by ocean model
          accum[i][j] = DEFAULT_ACCUMULATION_IN_OCEAN0;
          break;
        case 2:
          mask[i][j] = MASK_FLOATING;
          break;
        case 1:
        case 3:
          mask[i][j] = MASK_SHEET;
          // will determine below whether actually SHEET or DRAGGING
          break;
        default:
          SETERRQ(1,"invalid mask value in NetCDF initialization file");
      }
    }
  }
  // note ghosted values need to be communicated because of vote-by-neighbors
  // in IceModel::updateSurfaceElevationAndMask()
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  if (balVelRule != PETSC_TRUE) {
    ierr = PetscPrintf(grid.com, "done\n"); CHKERRQ(ierr);
    return 0;
  }
  
  // NOTE!!: from here on we are using mass balance velocities in the 
  //   computation of the mask:
  
  // compute deformational velocities (i.e. SIA and no Macayeal)
  const PetscTruth  saveUseMacayealVelocity = useMacayealVelocity;
  useMacayealVelocity = PETSC_FALSE;
  ierr = velocity(false); CHKERRQ(ierr);  // no need to update at depth; just
  // want ubar, vbar
  ierr = vertAveragedVelocityToRegular(); CHKERRQ(ierr); // communication here
  useMacayealVelocity = saveUseMacayealVelocity;
  
  // remove deformational from balance velocity to give sliding
  ierr = VecDuplicate(vh,&vsliding); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vsliding, &sliding); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbalvel, &balvel); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar deformationalMag =
        sqrt(ubar[i][j]*ubar[i][j] + vbar[i][j]*vbar[i][j]);
      sliding[i][j] = balvel[i][j] - deformationalMag;  // could be negative,
      // in which case get SHEET ...
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbalvel, &balvel); CHKERRQ(ierr);
  
  // now apply cutoff to determine mask
  const PetscScalar DEFAULT_MIN_SLIDING_FOR_MACAYEAL = 40.0; // m/a; FIXME: this is rigid
  const PetscScalar slideVelCutoff = DEFAULT_MIN_SLIDING_FOR_MACAYEAL / secpera;  

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // recall either Grounded or Rock became MASK_SHEET
      if (intMask(mask[i][j]) == MASK_SHEET) {
        if ((H[i][j] > 0.0) && (sliding[i][j] > slideVelCutoff)) {
          mask[i][j] = MASK_DRAGGING;
        } else {
          mask[i][j] = MASK_SHEET;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da2, vsliding, &sliding); CHKERRQ(ierr);
  ierr = VecDestroy(vsliding); CHKERRQ(ierr);

  // note ghosted values need to be communicated for next singleton removal
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  // now remove singleton and nearly singleton DRAGGING points, i.e. if 
  //   all BOX stencil neighbors are SHEET or if at most one is DRAGGING:
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) == MASK_DRAGGING) {
        const int neighmasksum = 
          modMask(mask[i-1][j+1]) + modMask(mask[i][j+1]) + modMask(mask[i+1][j+1]) +      
          modMask(mask[i-1][j])   +                       + modMask(mask[i+1][j])   +
          modMask(mask[i-1][j-1]) + modMask(mask[i][j-1]) + modMask(mask[i+1][j-1]);
        if (neighmasksum <= (7*MASK_SHEET + MASK_DRAGGING)) { 
          mask[i][j] = MASK_SHEET;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  // note ghosted values do need to be communicated because of vote-by-neighbors
  // in IceModel::updateSurfaceElevationAndMask()
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  ierr = PetscPrintf(grid.com, "done\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::bootstrapFromFile_netCDF_legacyAnt(const char *fname) {
  PetscErrorCode  ierr;
  int stat;
  char filename[PETSC_MAX_PATH_LEN];

  verbPrintf(2, grid.com, "Bootstrapping from Legacy Option\n");

  // The netCDF file has this physical extent
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  // FIXME: following is clearly tied to antarctica only!
  ierr = grid.rescale(280 * 20e3 / 2, 280 * 20e3 / 2, 5000); CHKERRQ(ierr);

  if (fname == NULL) { // This can be called from the driver for a default
    strcpy(filename, "pre_init.nc");
  } else {
    strcpy(filename, fname);
  }

  int ncid;
  int v_lon, v_lat, v_accum, v_h, v_H, v_bed, v_gl, v_T, v_ghf, v_uplift, v_balvel;
  if (grid.rank == 0) {
    stat = nc_open(filename, 0, &ncid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "lon", &v_lon); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "lat", &v_lat); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "ac", &v_accum); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "surf", &v_h); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "thk", &v_H); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "bed", &v_bed); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "gl", &v_gl); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "temps", &v_T); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "ghf", &v_ghf); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "uplift", &v_uplift); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "balvel", &v_balvel); CHKERRQ(nc_check(stat));
  }
  
  // COMMENT FIXME:  Next block of code uses a sequential Vec ("vzero") on processor 
  //   zero to take a NetCDF variable, which is a one-dimensional array representing 
  //   a two-dimensional (map-plane) quantity to the correct two-dimensional DA-based 
  //   Vec in IceModel.  A "scatter" to vector zero is needed because the entire 1D
  //   NetCDF variable is put into memory from the NetCDF variable on processor zero.
  //   Then on processor zero the VecSetValues puts the values in the global vec g2
  //   using the scatter context.  Then the values of g are put into local variables.
  //   All of the actual transfers, after the creation of the scatter context ctx, are
  //   done in ncVarToDAVec. 
  Vec vzero;
  VecScatter ctx;
  ierr = VecScatterCreateToZero(g2, &ctx, &vzero); CHKERRQ(ierr);  
  ierr = getIndZero(grid.da2, g2, vzero, ctx); CHKERRQ(ierr);
  // Now vzero contains indices to put processor zero quantities into
  // distributed global vectors.  These use g2 for scratch: the NetCDF variable
  // is put into an array on proc 0 then vzero is used to get the indices which
  // put the array into the scratch global Vec g2 and then the usual DA-based
  // global (g2) to local (vAccum, etc.) is done
  ierr = ncVarToDAVec(ncid, v_lon, grid.da2, vLongitude, g2, vzero);   CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_lat, grid.da2, vLatitude, g2, vzero);   CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_accum, grid.da2, vAccum, g2, vzero);   CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_h, grid.da2, vh, g2, vzero);       CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_H, grid.da2, vH, g2, vzero);       CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_bed, grid.da2, vbed, g2, vzero);     CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_gl, grid.da2, vMask, g2, vzero);    CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_T, grid.da2, vTs, g2, vzero);      CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_ghf, grid.da2, vGhf, g2, vzero);     CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_uplift, grid.da2, vuplift, g2, vzero);     CHKERRQ(ierr);

  // balvel is only locally used, so create and destroy here
  ierr = VecDuplicate(vh,&vbalvel); CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, v_balvel, grid.da2, vbalvel, g2, vzero);     CHKERRQ(ierr);
  
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }
  
  // At this point, the data still contains some missing data that we need to
  // fix.  Mostly, we will just fill in values.  Also some data has wrong units.
  // This is because the netCDF file should use intuitive values so that it can
  // be useful in its own right, however we use strictly SI units.  The method
  // below will fix this by scaling
  // Accumulation:    meters / a  -->  meters / s
  // bal vels:        meters / a  -->  meters / s
  // Temperature:     Celsius     -->  Kelvin
  // uplift:          meters / a  -->  meters / s
  ierr = cleanInputData_legacy(); CHKERRQ(ierr);

  // fill in temps at depth in reasonable way using surface temps and Ghf
  ierr = putTempAtDepth(); CHKERRQ(ierr);

#if 0
  PetscReal   maxH, maxh;
  ierr = VecMax(vH, PETSC_NULL, &maxH); CHKERRQ(ierr);
  ierr = VecMax(vh, PETSC_NULL, &maxh); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com, 
                      "properties of input data: Max(H) = %12.3f, Max(h) = %12.3f\n", 
                      maxH, maxh); CHKERRQ(ierr);

  PetscReal   maxMask, minMask;
  ierr = VecMax(vMask, PETSC_NULL, &maxMask); CHKERRQ(ierr);
  ierr = VecMin(vMask, PETSC_NULL, &minMask); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com, 
                      "properties of input data: Max(Mask) = %12.3f, Min(Mask) = %12.3f\n", 
                      maxMask, minMask); CHKERRQ(ierr);
#endif

  setConstantGrainSize(DEFAULT_GRAIN_SIZE);
  setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);

  ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);  
  // FIXME: vHmelt could probably be part of saved state.  In this case the best
  // procedure would be to *check* if vHmelt was saved, and if so to load it,
  // otherwise to set it to zero and report that.  Similar behavior appropriate
  // for many state variables.

  // for now: go ahead and create mask according to (balvel>cutoff) rule
  ierr = createMask_legacy(PETSC_TRUE); CHKERRQ(ierr);
  
  // vbalvel will not be used again ...
  ierr = VecDestroy(vbalvel); CHKERRQ(ierr);

  initialized_p = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceModel::reconfigure_legacy_Mbz() {
  PetscErrorCode  ierr;
  PetscTruth legacy_nc_Mbz;
  ierr = PetscOptionsHasName(PETSC_NULL, "-legacy_nc_Mbz", &legacy_nc_Mbz); CHKERRQ(ierr);
  if (legacy_nc_Mbz == PETSC_TRUE) {
    PetscInt    M, N, m, n;
    PetscScalar ***T, ***Tb, ***newTb;
    DA _da3b;
    Vec _vTb, _g3b;

    grid.p->Mbz++;
    ierr = DAGetInfo(grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                     PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
    ierr = DACreate3d(grid.com, DA_YZPERIODIC, DA_STENCIL_STAR, grid.p->Mbz, N, M,
                      1, n, m, 1, 1, PETSC_NULL, PETSC_NULL, PETSC_NULL, &_da3b); CHKERRQ(ierr);
    ierr = DACreateLocalVector(_da3b, &_vTb); CHKERRQ(ierr);
    ierr = DACreateGlobalVector(_da3b, &_g3b); CHKERRQ(ierr);

    ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
    ierr = DAVecGetArray(_da3b, _vTb, &newTb); CHKERRQ(ierr);
    for (PetscInt i = grid.xs; i < grid.xs + grid.xm; i++) {
      for (PetscInt j = grid.ys; j < grid.ys + grid.ym; j++) {
        for (PetscInt k = 0; k < grid.p->Mbz - 1; k++) {
          newTb[i][j][k] = Tb[i][j][k];
        }
        newTb[i][j][grid.p->Mbz - 1] = T[i][j][0];
      }
    }
    ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(_da3b, _vTb, &newTb); CHKERRQ(ierr);
    DADestroy(grid.da3b); grid.da3b = _da3b;
    VecDestroy(vTb); vTb = _vTb;
    VecDestroy(g3b); g3b = _g3b;
    // We never need the ghosted values of vTb, so we don't need LocalToLocal() here
  }
  return 0;
}
