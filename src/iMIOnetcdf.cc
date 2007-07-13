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
#include <cstdlib>
#include <petscda.h>
#include "iceModel.hh"

/* This first procedure does not actually refer to the NetCDF construct, but only to 
the scatter from processor zero.  It could, potentially, be used by other code which 
does the same scatter to processor zero. */
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


PetscErrorCode IceModel::putTempAtDepth() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **Ts, **Ghf, ***T, ***Tb;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar HH = H[i][j];
      const PetscInt    ks = static_cast<PetscInt>(floor(HH/grid.p->dz));
      const PetscScalar g = Ghf[i][j];
      // apply rule:  T(k) = Ts              above k=ks
      //              T(k) = Ts + alpha (H-z(k))^2 + beta (H-z(k))^4
      //                                  at k=0,1,...,ks
      //              Tb(l) = T(0) + (g/bedrock.k) * (Mbz - l) * dz
      //                                  at l=0,1,...,Mbz-1 
      // where alpha, beta are chosen so that dT/dz(z=0) = -g / ice.k,
      // and dT/dz(z=H/4) = -g / (2 ice.k)
      const PetscScalar beta = (4.0/21.0)
        * (g / (2.0 * ice.k * HH * HH * HH));
      const PetscScalar alpha = (g / (2.0 * HH * ice.k)) - 2.0 * HH * HH * beta;
      for (PetscInt k=ks; k<grid.p->Mz; k++) {
        T[i][j][k] = Ts[i][j];
      }
      for (PetscInt k=0; k<ks; k++) {
        const PetscScalar depth = HH - static_cast<PetscScalar>(k) * grid.p->dz;
        const PetscScalar depth2 = depth*depth;
        T[i][j][k] = Ts[i][j] + alpha * depth2 + beta * depth2 * depth2;
      }
      for (PetscInt kb=0; kb<grid.p->Mbz; kb++)
        Tb[i][j][kb] = T[i][j][0] + (Ghf[i][j]/bedrock.k) * 
          static_cast<PetscScalar>(grid.p->Mbz - kb - 1) * grid.p->dz;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3b, vTb, INSERT_VALUES, vTb); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3b, vTb, INSERT_VALUES, vTb); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::setAccumInOcean() {
  PetscErrorCode ierr;
  PetscScalar **mask, **accum;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask[i][j] == MASK_FLOATING || mask[i][j] == MASK_FLOATING_OCEAN0) {
        accum[i][j] = DEFAULT_ACCUMULATION_IN_OCEAN0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  // no communication needed for accum
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

#include <netcdf.h>
#include "nc_util.hh"

PetscErrorCode IceModel::nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}

PetscErrorCode IceModel::bootstrapFromFile_netCDF(const char *fname) {
  PetscErrorCode  ierr;
  int stat;
  PetscReal lonExists, latExists, accumExists, hExists, HExists, bExists,
      TsExists, ghfExists, upliftExists, balvelExists,
      xExists, yExists, HmeltExists, tExists;
  PetscReal glonExists=0, glatExists=0, gaccumExists=0, ghExists=0, gHExists=0, gbExists=0,
      gTsExists=0, gghfExists=0, gupliftExists=0, gbalvelExists=0,
      gxExists=0, gyExists=0, gHmeltExists=0, gtExists=0;

  verbPrintf(2, grid.com, "bootstrapping by PISM default method from file %s\n",fname);

  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);

  if (fname == NULL) {
    SETERRQ(1, "No file name given for bootstrapping\n");
  }

  // determine if variables exist in bootstrapping file
  int ncid;
  int v_lon, v_lat, v_accum, v_h, v_H, v_bed, v_Ts, v_ghf, v_uplift,
      v_balvel, v_x, v_y, v_Hmelt, v_t;
  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "x", &v_x); gxExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "y", &v_y); gyExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "lon", &v_lon); glonExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "lat", &v_lat); glatExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "accum", &v_accum); gaccumExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "h", &v_h); ghExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "H", &v_H); gHExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "b", &v_bed); gbExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "Ts", &v_Ts); gTsExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "ghf", &v_ghf); gghfExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "uplift", &v_uplift); gupliftExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "balvel", &v_balvel); gbalvelExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "Hmelt", &v_Hmelt); gHmeltExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "t", &v_t); gtExists = stat == NC_NOERR;
  }
  ierr = PetscGlobalMax(&gxExists, &xExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gyExists, &yExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&glonExists, &lonExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&glatExists, &latExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gaccumExists, &accumExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&ghExists, &hExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gHExists, &HExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gbExists, &bExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gTsExists, &TsExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gghfExists, &ghfExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gupliftExists, &upliftExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gbalvelExists, &balvelExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gHmeltExists, &HmeltExists, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&gtExists, &tExists, grid.com); CHKERRQ(ierr);
  

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
  PetscScalar first, last;
  PetscScalar x_scale = 750.0e3, y_scale = 750.0e3, z_scale = 5000.0;
  if (xExists) {
    ierr = getFirstLast(ncid, v_x, &first, &last);
    x_scale = (last - first) / 2.0;
    // verbPrintf(1,grid.com,"   for x: last,first,x_scale = %f,%f,%f\n",last,first,x_scale);
  } else verbPrintf(2, grid.com, 
         "  WARNING: variable x(x) not found in file. Proceeding with default interval [-Lx,Lx] = [-%f,%f]\n",
         x_scale,x_scale);
  if (yExists) {
    ierr = getFirstLast(ncid, v_y, &first, &last);
    y_scale = (last - first) / 2.0;
    // verbPrintf(1,grid.com,"   for y: last,first,y_scale = %f,%f,%f\n",last,first,y_scale);
  } else verbPrintf(2, grid.com, 
         "  WARNING: variable y(y) not found in file. Proceeding with defaults [-Ly,Ly] = [-%f,%f]\n",
         y_scale,y_scale);
  verbPrintf(2, grid.com, 
         "  using default value Lz=%f for vertical extent of computational box for ice\n", z_scale);
  ierr = grid.rescale(x_scale, y_scale, z_scale); CHKERRQ(ierr);

  if (lonExists) {
    ierr = ncVarToDAVec(ncid, v_lon, grid.da2, vLongitude, g2, vzero);
    CHKERRQ(ierr);
  } else verbPrintf(3, grid.com, "  continuing without lon\n"); 
  if (latExists) {
    ierr = ncVarToDAVec(ncid, v_lat, grid.da2, vLatitude, g2, vzero);
    CHKERRQ(ierr);
  } else verbPrintf(3, grid.com, "  continuing without lat\n"); 
  if (accumExists) {
    ierr = ncVarToDAVec(ncid, v_accum, grid.da2, vAccum, g2, vzero);
    CHKERRQ(ierr);
  } else {
    verbPrintf(2, grid.com, "  accum not found. Filling in with default value: %f\n", DEFAULT_ACCUM_VALUE_MISSING); 
    ierr = VecSet(vAccum, DEFAULT_ACCUM_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (hExists)  
    verbPrintf(2, grid.com, "  WARNING: ignoring values found for surface elevation h and using h = b + H\n");
  if (HExists) {
    ierr = ncVarToDAVec(ncid, v_H, grid.da2, vH, g2, vzero);
    CHKERRQ(ierr);
  } else {
    verbPrintf(2, grid.com, "  WARNING: thickness H not found. Filling in with default value: %f\n",
               DEFAULT_H_VALUE_MISSING);
    ierr = VecSet(vH, DEFAULT_H_VALUE_MISSING); CHKERRQ(ierr); 
  }
  if (bExists) {
    ierr = ncVarToDAVec(ncid, v_bed, grid.da2, vbed, g2, vzero);
    CHKERRQ(ierr);
  } else {
    verbPrintf(2, grid.com, "  WARNING: bedrock elevation b not found. Filling in with default value: %f\n",
               DEFAULT_BED_VALUE_MISSING); 
    ierr = VecSet(vbed, DEFAULT_BED_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (TsExists) {
    ierr = ncVarToDAVec(ncid, v_Ts, grid.da2, vTs, g2, vzero); CHKERRQ(ierr);
  } else {
    verbPrintf(2, grid.com, "  WARNING: surface temperature Ts not found. Filling in with default value: %f\n",
               DEFAULT_SURF_TEMP_VALUE_MISSING);
    ierr = VecSet(vTs, DEFAULT_SURF_TEMP_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (ghfExists) {
    ierr = ncVarToDAVec(ncid, v_ghf, grid.da2, vGhf, g2, vzero); CHKERRQ(ierr);
  } else {
    verbPrintf(2, grid.com, "  WARNING: geothermal flux ghf not found. Filling in with default value: %f\n",
               DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING); 
    ierr = VecSet(vGhf, DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (upliftExists) {
    ierr = ncVarToDAVec(ncid, v_uplift, grid.da2, vuplift, g2, vzero); 
    CHKERRQ(ierr);
  } else {
    verbPrintf(3, grid.com, "  uplift not found. Filling with zero\n"); 
    ierr = VecSet(vuplift, 0.0); CHKERRQ(ierr);
  }
  verbPrintf(2, grid.com, "  determining mask by floating criterion; grounded ice marked as SIA (=1)\n");
  ierr = setMaskSurfaceElevation_bootstrap(); CHKERRQ(ierr);
  verbPrintf(3, grid.com, 
             "  setting accumulation in ice shelf-free ocean to default value %f\n",
             DEFAULT_ACCUMULATION_IN_OCEAN0 * secpera);
  ierr = setAccumInOcean(); CHKERRQ(ierr);
  if (HmeltExists) {
    ierr = ncVarToDAVec(ncid, v_Hmelt, grid.da2, vHmelt, g2, vzero); CHKERRQ(ierr);
  } else {
    verbPrintf(2, grid.com, "  Hmelt not found. Filling with zero\n");
    ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);  
  }
  if (tExists) {
    const size_t idx = 0;
    float begin_t;
    nc_get_var1_float(ncid, v_t, &idx, &begin_t);
    setInitialAgeYears(begin_t);
  } else {
    verbPrintf(2, grid.com, "  time t not found. Using t=%f years\n",
               DEFAULT_INITIAL_AGE_YEARS);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
  }
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }

  // fill in temps at depth in reasonable way using surface temps and Ghf
  ierr = putTempAtDepth(); CHKERRQ(ierr);

  setConstantGrainSize(DEFAULT_GRAIN_SIZE);

  verbPrintf(2, grid.com, "bootstrapping done\n");
  initialized_p = PETSC_TRUE;
  return 0;
}

// This function is only here until the new and better bootstrap function is
// complete.
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

// change to use MPIBcast
PetscErrorCode IceModel::ncVarBcastVec(int ncid, int vid, Vec *vecg) {
  PetscErrorCode ierr;
  int stat;
  size_t M;
  float *f = NULL;

  if (grid.rank == 0) {
    int dimids[NC_MAX_VAR_DIMS];
    int ndims, natts;
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 1) {
      SETERRQ2(1, "ncVarBcastDaVec: number of dimensions = %d for %s\n",
               ndims, name);
    }
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));
    f = new float[M];
    stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
  }

  // broadcast the length
  ierr = MPI_Bcast(&M, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
 
  // if you're not rank 0, you still need to create the array
  if (grid.rank != 0){
    f = new float[M];
  }
  ierr = MPI_Bcast(f, M, MPI_FLOAT, 0, grid.com); CHKERRQ(ierr);

  VecCreateSeq(grid.com, M, vecg);

  for (int x=0; x<(int)M; x++) {
    ierr = VecSetValue(*vecg, x, f[x], INSERT_VALUES); CHKERRQ(ierr);
  }
  
  ierr = VecAssemblyBegin(*vecg); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*vecg); CHKERRQ(ierr)

  delete [] f;

  return 0;
}

PetscErrorCode IceModel::ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                                      Vec vecg, Vec vindzero) {
  PetscErrorCode  ierr;
  int stat;
  PetscScalar **ind;
  float       *f = NULL;
  int         *g = NULL;

  if (grid.rank == 0) {
    int dimids[NC_MAX_VAR_DIMS];
    int ndims, natts;
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 2) {
      SETERRQ2(1, "ncVarToDaVec: number of dimensions = %d for %s\n",
               ndims, name);
    }

    // In the netCDF file,
    // we index 0:M in the x direction and 0:N in the y direction.  Such a
    // location $(i,j) \in [0,M] \times [0,N]$ is addressed as [i*N + j]
    size_t M, N;
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));
    stat = nc_inq_dimlen(ncid, dimids[1], &N); CHKERRQ(nc_check(stat));
    
    switch (xtype) {
      case NC_INT:
        g = new int[M*N];
        stat = nc_get_var_int(ncid, vid, g); CHKERRQ(nc_check(stat));
        break;
      case NC_FLOAT:
        f = new float[M*N];
        stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
        break;
      default:
        SETERRQ1(1, "NC_VAR `%s' not of type NC_INT or NC_FLOAT.\n", name);
    }

    ierr = VecGetArray2d(vindzero, grid.p->Mx, grid.p->My, 0, 0, &ind); CHKERRQ(ierr);
    // ierr = VecGetArray(vzero, &w_flat); CHKERRQ(ierr);
    
    // netCDF concepts of $\Delta x$ and $\Delta y$
    // We have rescaled the grid early in initFromFile_netCDF()
    // to match the physical extent of the netCDF file.
    const float ncdx = 2 * grid.p->Lx / (M - 1);
    const float ncdy = 2 * grid.p->Ly / (N - 1);

    for (PetscInt i=0; i < grid.p->Mx; i++) {
      for (PetscInt j=0; j < grid.p->My; j++) {
        const float x = grid.p->dx * (i - grid.p->Mx/2);
        const float y = grid.p->dy * (j - grid.p->My/2);
        if ((PetscAbs(x) > grid.p->Lx) || (PetscAbs(y) > grid.p->Ly)) {
          SETERRQ1(1, "ncVarToDaVec: x=%f not in bounds.  Grid corrupted.\n", x);
        }
            
        const float ii = M / 2 + x / ncdx;
        const float jj = N / 2 + y / ncdy;
        // These live in [0,1)
        const float xx = ii - floor(ii);
        const float yy = jj - floor(jj);
        // Define weights for bilinear interpolation
        const float w11 = (1 - xx) * (1 - yy);
        const float w12 = xx * (1 - yy);
        const float w21 = (1 - xx) * yy;
        const float w22 = xx * yy;
        // Locations to sample from.  These should be on the grid.
        const int i1 = int(floor(ii)) % M;
        const int i2 = int(ceil(ii)) % M;
        const int j1 = int(floor(jj)) % N;
        const int j2 = int(ceil(jj)) % N;

        PetscScalar val;
        if (g != NULL) { // We have an integer array
          val = (w11 * g[i1*N + j1] + w12 * g[i1*N + j2]
                 + w21 * g[i2*N + j1] + w22 * g[i2*N + j2]);
        } else if (f != NULL) { // It is a float array
          val = (w11 * f[i1*N + j1] + w12 * f[i1*N + j2]
                 + w21 * f[i2*N + j1] + w22 * f[i2*N + j2]);
        } else {
          SETERRQ(1, "This should not happen");
        }
        
        // The backward indexing is merely to make the plots look upright with
        // the default plotting methods.  When I can make the axes work in a
        // sane manner, this can be improved.
        ierr = VecSetValue(vecg, (PetscInt) ind[grid.p->Mx - 1 - i][j],
                           val, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
   
    if (f != NULL) delete [] f;
    if (g != NULL) delete [] g;
    
    ierr = VecRestoreArray2d(vindzero, grid.p->Mx, grid.p->My, 0, 0, &ind); CHKERRQ(ierr);
  }

  ierr = VecAssemblyBegin(vecg); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vecg); CHKERRQ(ierr);

  ierr = DAGlobalToLocalBegin(da, vecg, INSERT_VALUES, vecl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, vecg, INSERT_VALUES, vecl); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::getFirstLast(int ncid, int vid, PetscScalar *gfirst, PetscScalar *glast) {
  int stat;
  PetscScalar first, last;
  float       *f = NULL;
  int         *g = NULL;

  if (grid.rank == 0) {
    int dimids[NC_MAX_VAR_DIMS];
    int ndims, natts;
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 1) {
      SETERRQ2(1, "getFirstLast: number of dimensions = %d for %s\n",
               ndims, name);
    }

    // In the netCDF file,
    // we index 0:M in the x direction and 0:N in the y direction.  Such a
    // location $(i,j) \in [0,M] \times [0,N]$ is addressed as [i*N + j]
    size_t M;
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));

    switch (xtype) {
      case NC_INT:
        g = new int[M];
        stat = nc_get_var_int(ncid, vid, g); CHKERRQ(nc_check(stat));
        break;
      case NC_FLOAT:
        f = new float[M];
        stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
        break;
      default:
        SETERRQ1(1, "NC_VAR `%s' not of type NC_INT or NC_FLOAT.\n", name);
    }

    if (g != NULL) {
      first = g[0];
      last = g[M-1];
    } else if (f != NULL) {
      first = f[0];
      last = f[M-1];
    } else {
      SETERRQ(1, "This should not happen.\n");
    }
  } else {
    first = 1.0e30;
    last = -1.0e30;
  }
  
  PetscGlobalMin(&first,gfirst,grid.com);  
  PetscGlobalMax(&last,glast,grid.com);
  return 0;
}

PetscErrorCode IceModel::setMaskSurfaceElevation_bootstrap() {
    PetscErrorCode ierr;
  PetscScalar **h, **bed, **H, **mask;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // take this opportunity to check that H[i][j] >= 0
      if (H[i][j] < 0.0) {
        SETERRQ2(1,"Thickness negative at point i=%d, j=%d",i,j);
      }
      
      if (H[i][j] < 0.001) {  // if no ice
        if (bed[i][j] < 0.0) {
          h[i][j] = 0.0;
          mask[i][j] = MASK_FLOATING_OCEAN0;
        } else {
          h[i][j] = bed[i][j];
          mask[i][j] = MASK_SHEET;
        } 
      } else { // if some ice thickness then check floating criterion
        const PetscScalar hgrounded = bed[i][j] + H[i][j];
        const PetscScalar hfloating = (1-ice.rho/ocean.rho) * H[i][j];
        // check whether you are actually floating or grounded
        if (hgrounded > hfloating) {
          h[i][j] = hgrounded; // actually grounded so set h
          mask[i][j] = MASK_SHEET;
        } else {
          h[i][j] = hfloating; // actually floating so update h
          mask[i][j] = MASK_FLOATING;
        }
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  // go ahead and communicate mask and surface elev now; may be redundant communication?
  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::dumpToFile_netCDF(const char *fname) {
  PetscErrorCode ierr;

#include "netcdf/write_attributes.c"

  if (grid.rank == 0) {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "history",
                           strlen(grid.p->history), grid.p->history);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
    
    stat = nc_enddef (ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    // float t = grid.p->year * secpera;
    double t = grid.p-> year * secpera;
    size_t zero = 0;
    // stat = nc_put_var1_float(ncid, t_id, &zero, &t); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_var1_double(ncid, t_id, &zero, &t); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    ierr = put_dimension_regular(ncid, x_id, x_len, -grid.p->Lx, grid.p->dx); CHKERRQ(ierr);
    ierr = put_dimension_regular(ncid, y_id, y_len, -grid.p->Ly, grid.p->dy); CHKERRQ(ierr);
    ierr = put_dimension_regular(ncid, z_id, z_len, 0, grid.p->dz); CHKERRQ(ierr);
    ierr = put_dimension_regular(ncid, zb_id, zb_len, -grid.p->Lbz, grid.p->dz); CHKERRQ(ierr);
  }

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.p->Mbz}; // Count local block: bed

  // Allocate some memory.  We will assume that vectors based on grid.da3 are the largest.
  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

  // these are treated like 2-D constant quantities
  ierr = put_local_var(&grid, ncid, lon_id, NC_FLOAT, grid.da2, vLongitude, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, lat_id, NC_FLOAT, grid.da2, vLatitude, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D model quantities
  ierr = put_local_var(&grid, ncid, mask_id, NC_BYTE, grid.da2, vMask, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, h_id, NC_FLOAT, grid.da2, vh, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, H_id, NC_FLOAT, grid.da2, vH, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, Hmelt_id, NC_FLOAT, grid.da2, vHmelt, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, b_id, NC_FLOAT, grid.da2, vbed, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, dbdt_id, NC_FLOAT, grid.da2, vuplift, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D climate quantities
  ierr = put_local_var(&grid, ncid, Ts_id, NC_FLOAT, grid.da2, vTs, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, ghf_id, NC_FLOAT, grid.da2, vGhf, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, accum_id, NC_FLOAT, grid.da2, vAccum, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
/*                       
  // 2-D constant quantities
  ierr = put_local_var(&grid, ncid, Ts_id, NC_FLOAT, grid.da2, vTs, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, ghf_id, NC_FLOAT, grid.da2, vGhf, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, accum_id, NC_FLOAT, grid.da2, vAccum, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
*/
  // 3-D model quantities
  ierr = put_local_var(&grid, ncid, T_id, NC_FLOAT, grid.da3, vT, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, Tb_id, NC_FLOAT, grid.da3b, vTb, g3b,
                       s, cb, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, age_id, NC_FLOAT, grid.da3, vtau, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);

  // We are done with these buffers
  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}


PetscErrorCode IceModel::initFromFile_netCDF(const char *fname) {
  PetscErrorCode  ierr;
  size_t      dim[5];
  float       bdy[7];
  double 	  bdy_time;
  int         ncid, stat;
  PetscInt    ignor;
  PetscTruth  M_Set;

  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  ierr = get_dimensions(ncid, dim, bdy, &bdy_time, grid.com); CHKERRQ(ierr);

  // grid.p->year = bdy[0] / secpera;
  grid.p->year = bdy_time / secpera;
  grid.p->Mx = dim[1];
  grid.p->My = dim[2];
  grid.p->Mz = dim[3];
  grid.p->Mbz = dim[4];
  // options read simply to warn user that they have been ignored
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mx", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,"WARNING: user option -Mx ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-My", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,"WARNING: user option -My ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mz", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,"WARNING: user option -Mz ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mbz", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,"WARNING: user option -Mbz ignored; value read from file %s\n", fname); CHKERRQ(ierr);
  }
  
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);

  ierr = grid.rescale(-bdy[1], -bdy[3], bdy[6]); CHKERRQ(ierr);
  // These lines remain only to confirm the use of the fields.  If we have a
  // more flexible gridding capability, the information in the unused fields may
  // be useful.
  //   grid.p->Lx = -bdy[1]; // == bdy[2]
  //   grid.p->Ly = -bdy[3]; // == bdy[4]
  //   grid.p->Lbz = -bdy[5];
  //   grid.p->Lz = bdy[6];
  // note user setting of -Lx,-Ly,-Lz will overwrite these settings from file

  // set IceModel::startYear, IceModel::endYear, grid.p->year, but respecting grid.p->year
  // which came from -if file, _unless_ -ys set by user
  ierr = setStartRunEndYearsFromOptions(PETSC_TRUE);  CHKERRQ(ierr);

  // Time to compute what we need.
  int s[] = {dim[0] - 1, grid.xs, grid.ys, 0};   // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.p->Mbz}; // Count local block: bed

  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

  // these are treated like 2-D constant quantities
  ierr = get_local_var(&grid, ncid, "lon", NC_FLOAT, grid.da2, vLongitude, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "lat", NC_FLOAT, grid.da2, vLatitude, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D model quantities
  ierr = get_local_var(&grid, ncid, "mask", NC_BYTE, grid.da2, vMask, g2, s, c, 3,
                       a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "h", NC_FLOAT, grid.da2, vh, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "H", NC_FLOAT, grid.da2, vH, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "Hmelt", NC_FLOAT, grid.da2, vHmelt, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "b", NC_FLOAT, grid.da2, vbed, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "dbdt", NC_FLOAT, grid.da2, vuplift, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  // 2-D climate quantities
  ierr = get_local_var(&grid, ncid, "Ts", NC_FLOAT, grid.da2, vTs, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "ghf", NC_FLOAT, grid.da2, vGhf, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "accum", NC_FLOAT, grid.da2, vAccum, g2,
                       s, c, 3, a_mpi, max_a_len); CHKERRQ(ierr);
/*
  // 2-D constant quantities
  ierr = get_local_var(&grid, ncid, "Ts", NC_FLOAT, grid.da2, vTs, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "ghf", NC_FLOAT, grid.da2, vGhf, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "accum", NC_FLOAT, grid.da2, vAccum, g2,
                       &s[1], &c[1], 2, a_mpi, max_a_len); CHKERRQ(ierr);
*/
  // 3-D model quantities
  ierr = get_local_var(&grid, ncid, "T", NC_FLOAT, grid.da3, vT, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "Tb", NC_FLOAT, grid.da3b, vTb, g3b,
                       s, cb, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "age", NC_FLOAT, grid.da3, vtau, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);

  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

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

  initialized_p = PETSC_TRUE;
  return 0;
}

