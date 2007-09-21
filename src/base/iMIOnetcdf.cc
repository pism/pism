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


PetscErrorCode IceModel::putTempAtDepth() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **b, **Ts, **Ghf, ***T, ***Tb;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar HH = H[i][j];
      const PetscInt    ks = static_cast<PetscInt>(floor(HH/grid.p->dz));
      if (ks >= grid.p->Mz) {
        ierr = verbPrintf(1,grid.com,
                 "[[error LOCATION: i, j, ks, H = %5d %5d %5d %10.2f]]\n",
                 i, j, ks, H[i][j]);  CHKERRQ(ierr); 
        SETERRQ(1,"Vertical grid exceeded");
      }
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
        const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * depth;
        const PetscScalar d2 = depth*depth;
        T[i][j][k] = PetscMin(Tpmp,Ts[i][j] + alpha * d2 + beta * d2 * d2);
      }
      PetscScalar T_top_bed = T[i][j][0];
      // if floating then top of bedrock sees ocean
      const PetscScalar floating_base = - (ice.rho/ocean.rho) * H[i][j];
      if (b[i][j] < floating_base - 1.0)
        T_top_bed = ice.meltingTemp;
      for (PetscInt kb=0; kb<grid.p->Mbz; kb++) {
        Tb[i][j][kb] = T_top_bed + (Ghf[i][j]/bedrock.k) * 
          static_cast<PetscScalar>(grid.p->Mbz - kb - 1) * grid.p->dz;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
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
      if (mask[i][j] == MASK_FLOATING_OCEAN0) {
        accum[i][j] = DEFAULT_ACCUMULATION_IN_OCEAN0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  // no communication needed for accum
  return 0;
}


PetscErrorCode IceModel::nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}


PetscErrorCode IceModel::bootstrapFromFile_netCDF(const char *fname) {
  PetscErrorCode  ierr;
  int stat;
  PetscInt psExists=0, lonExists=0, latExists=0, accumExists=0, hExists=0, HExists=0, bExists=0,
      TsExists=0, ghfExists=0, upliftExists=0, balvelExists=0,
      xExists=0, yExists=0, HmeltExists=0, tExists=0;

  ierr = verbPrintf(2, grid.com, "bootstrapping by PISM default method from file %s\n",fname); 
           CHKERRQ(ierr);
  if (fname == NULL) {
    SETERRQ(1, "No file name given for bootstrapping\n");
  }

  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);

  // determine if variables exist in bootstrapping file
  int ncid;
  int v_ps, v_lon, v_lat, v_accum, v_h, v_H, v_bed, v_Ts, v_ghf, v_uplift,
      v_balvel, v_x, v_y, v_Hmelt, v_t;
  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "polar_stereographic", &v_ps); psExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "x", &v_x); xExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "y", &v_y); yExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "t", &v_t); tExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "lon", &v_lon); lonExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "lat", &v_lat); latExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "accum", &v_accum); accumExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "h", &v_h); hExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "H", &v_H); HExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "b", &v_bed); bExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "Ts", &v_Ts); TsExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "ghf", &v_ghf); ghfExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "uplift", &v_uplift); upliftExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "balvel", &v_balvel); balvelExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "Hmelt", &v_Hmelt); HmeltExists = stat == NC_NOERR;
  }
  // broadcast the existence flags
  ierr = MPI_Bcast(&psExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&xExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&yExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&tExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&lonExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&latExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&accumExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&hExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&HExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&bExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&TsExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ghfExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&upliftExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&balvelExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&HmeltExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  
  // if polar_stereographic variable exists, then store attributes for writing out 
  // at end of run
  if (psExists) {
    if (grid.rank == 0) {
      stat = nc_get_att_double(ncid, v_ps, "straight_vertical_longitude_from_pole",
                              &psParams.svlfp); CHKERRQ(nc_check(stat));
      stat = nc_get_att_double(ncid, v_ps, "latitude_of_projection_origin",
                              &psParams.lopo); CHKERRQ(nc_check(stat));
      stat = nc_get_att_double(ncid, v_ps, "standard_parallel",
                              &psParams.sp); CHKERRQ(nc_check(stat));
    }
    ierr = MPI_Bcast(&psParams.svlfp, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&psParams.lopo, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&psParams.sp, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    ierr = verbPrintf(3,grid.com,
               "  polar stereographic found: svlfp = %6.2f, lopo = %6.2f, sp = %6.2f\n",
               psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr); 
  } else {
    ierr = verbPrintf(3,grid.com,
               "  polar stereo not found, using defaults: svlfp=%6.2f, lopo=%6.2f, sp=%6.2f\n",
               psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr); 
  }

  // set time (i.e. grid.p->year, IceModel::startYear, and IceModel::endYear)
  if (tExists) {
    const size_t idx = 0;
    float begin_t;
    if (grid.rank == 0) {
      nc_get_var1_float(ncid, v_t, &idx, &begin_t);
    }
    ierr = MPI_Bcast(&begin_t, 1, MPI_FLOAT, 0, grid.com); CHKERRQ(ierr);
    grid.p->year = begin_t / secpera;
    ierr = verbPrintf(3, grid.com, 
            "  time t found in bootstrap file; using to set current year\n"); CHKERRQ(ierr);
    ierr = setStartRunEndYearsFromOptions(PETSC_TRUE); CHKERRQ(ierr);
  } else {
    // set the current year from options or defaults only
    ierr = setStartRunEndYearsFromOptions(PETSC_FALSE); CHKERRQ(ierr);
  }

  // set grid x,y,z from bootstrap file
  PetscScalar first, last;
  PetscScalar x_scale = 750.0e3, y_scale = 750.0e3, z_scale = 4000.0; // just defaults
  if (xExists) {
    ierr = getFirstLast(ncid, v_x, &first, &last); CHKERRQ(ierr);
    x_scale = (last - first) / 2.0;
  } else {
    ierr = verbPrintf(2, grid.com,"  WARNING: variable x(x) not found in bootstrap file."
         " Proceeding with default interval [-Lx,Lx] = [-%f,%f]\n",x_scale,x_scale); CHKERRQ(ierr);
  }
  if (yExists) {
    ierr = getFirstLast(ncid, v_y, &first, &last); CHKERRQ(ierr);
    y_scale = (last - first) / 2.0;
  } else {
    ierr = verbPrintf(2, grid.com,"  WARNING: variable y(y) not found in bootstrap file."
         " Proceeding with default interval [-Ly,Ly] = [-%f,%f]\n",y_scale,y_scale); CHKERRQ(ierr);
  }
  ierr = grid.rescale(x_scale, y_scale, z_scale); CHKERRQ(ierr);
  // runtime options take precedence in setting of -Lx,-Ly,-Lz *including*
  // if initialization is from an input file
  PetscTruth LxSet, LySet, LzSet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lx", &x_scale, &LxSet); CHKERRQ(ierr);
  if (LxSet == PETSC_TRUE) {
    ierr = grid.rescale(x_scale*1000.0, grid.p->Ly, grid.p->Lz); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Ly", &y_scale, &LySet); CHKERRQ(ierr);
  if (LySet == PETSC_TRUE) {
    ierr = grid.rescale(grid.p->Lx, y_scale*1000.0, grid.p->Lz); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lz", &z_scale, &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_TRUE) {
    ierr = grid.rescale(grid.p->Lx, grid.p->Ly, z_scale); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, 
         "  using default value Lz=%f for vertical extent of computational box for ice\n",
         z_scale); CHKERRQ(ierr);
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
  if (lonExists) {
    ierr = ncVarToDAVec(ncid, v_lon, grid.da2, vLongitude, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  continuing without lon\n");  CHKERRQ(ierr);
  }
  if (latExists) {
    ierr = ncVarToDAVec(ncid, v_lat, grid.da2, vLatitude, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  continuing without lat\n");  CHKERRQ(ierr);
  }
  if (accumExists) {
    ierr = ncVarToDAVec(ncid, v_accum, grid.da2, vAccum, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
               "  accum not found; using default %7.2f m/a\n",
               DEFAULT_ACCUM_VALUE_MISSING * secpera);  CHKERRQ(ierr);
    ierr = VecSet(vAccum, DEFAULT_ACCUM_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (HExists) {
    ierr = ncVarToDAVec(ncid, v_H, grid.da2, vH, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
               "  WARNING: thickness H not found; using default %8.2f\n",
               DEFAULT_H_VALUE_MISSING); CHKERRQ(ierr);
    ierr = VecSet(vH, DEFAULT_H_VALUE_MISSING); CHKERRQ(ierr); 
  }
  if (bExists) {
    ierr = ncVarToDAVec(ncid, v_bed, grid.da2, vbed, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
             "  WARNING: bedrock elevation b not found; using default %f\n",
             DEFAULT_BED_VALUE_MISSING);  CHKERRQ(ierr);
    ierr = VecSet(vbed, DEFAULT_BED_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (hExists) {
    ierr = verbPrintf(2, grid.com, 
          "  WARNING: ignoring values found for surface elevation h and using h = b + H\n");
          CHKERRQ(ierr);
  }
  if (TsExists) {
    ierr = ncVarToDAVec(ncid, v_Ts, grid.da2, vTs, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
             "  WARNING: surface temperature Ts not found; using default %7.2f K\n",
             DEFAULT_SURF_TEMP_VALUE_MISSING); CHKERRQ(ierr);
    ierr = VecSet(vTs, DEFAULT_SURF_TEMP_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (ghfExists) {
    ierr = ncVarToDAVec(ncid, v_ghf, grid.da2, vGhf, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, 
             "  WARNING: geothermal flux ghf not found; using default %6.3f W/m^2\n",
             DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING);  CHKERRQ(ierr);
    ierr = VecSet(vGhf, DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (upliftExists) {
    ierr = ncVarToDAVec(ncid, v_uplift, grid.da2, vuplift, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  uplift not found. Filling with zero\n");  CHKERRQ(ierr);
    ierr = VecSet(vuplift, 0.0); CHKERRQ(ierr);
  }
  if (HmeltExists) {
    ierr = ncVarToDAVec(ncid, v_Hmelt, grid.da2, vHmelt, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  Hmelt not found. Filling with zero\n"); CHKERRQ(ierr);
    ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);  
  }
  // done with scatter context
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);
  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }
  ierr = verbPrintf(3, grid.com, "  done reading .nc file\n"); CHKERRQ(ierr);

  ierr = cleanJustNan_legacy(); CHKERRQ(ierr);  // needed for ant_init.nc; should be harmless for other

  ierr = verbPrintf(3, grid.com, 
            "  determining mask by floating criterion; grounded ice marked as SIA (=1)\n");
            CHKERRQ(ierr);
  ierr = setMaskSurfaceElevation_bootstrap(); CHKERRQ(ierr);
  ierr = verbPrintf(3, grid.com, 
             "  setting accumulation in ice shelf-free ocean to default value %6.2f m/a\n",
             DEFAULT_ACCUMULATION_IN_OCEAN0 * secpera); CHKERRQ(ierr);
  ierr = setAccumInOcean(); CHKERRQ(ierr);
  
  // fill in temps at depth in reasonable way using surface temps and Ghf
  ierr = verbPrintf(2, grid.com, 
             "  filling in temperatures at depth using surface temperatures and quartic guess\n"); CHKERRQ(ierr);
  ierr = putTempAtDepth(); CHKERRQ(ierr);
  // fill in other 3D fields
  setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
  setConstantGrainSize(DEFAULT_GRAIN_SIZE);

  verbPrintf(2, grid.com, "bootstrapping done\n");
  initialized_p = PETSC_TRUE;
  return 0;
}

PetscErrorCode IceModel::ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                                      Vec vecg, Vec vindzero) {
  PetscErrorCode  ierr;
  MaskInterp defaultmasktool;
  
  defaultmasktool.number_allowed = 4;
  defaultmasktool.allowed_levels[0] = MASK_SHEET;
  defaultmasktool.allowed_levels[1] = MASK_DRAGGING;
  defaultmasktool.allowed_levels[2] = MASK_FLOATING;
  defaultmasktool.allowed_levels[3] = MASK_FLOATING_OCEAN0;
  ierr = ncVarToDAVec(ncid,vid,da,vecl,vecg,vindzero,defaultmasktool); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                                      Vec vecg, Vec vindzero, MaskInterp masktool) {
  PetscErrorCode  ierr;
  int stat;
  PetscScalar **ind;
  float       *f = NULL;
  int         *g = NULL;

  if (masktool.number_allowed < 1) {
    SETERRQ(99, "ncvarToDAVec: number of allowed levels in masktool must be at least one");
  }
  if (grid.rank == 0) {
    int ndims, natts, dimids[NC_MAX_VAR_DIMS];
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 2) {
      SETERRQ2(1, "ncVarToDaVec: number of dimensions = %d for %s\n",ndims, name);
    }
    // in netCDF file, we index 0:M in the x direction and 0:N in the y direction
    // location (i,j) is addressed as [i*N + j]
    size_t M, N;
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));
    stat = nc_inq_dimlen(ncid, dimids[1], &N); CHKERRQ(nc_check(stat));
    
    switch (xtype) {
      case NC_FLOAT:
        f = new float[M*N];
        stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
        break;
      case NC_INT:
        g = new int[M*N];
        stat = nc_get_var_int(ncid, vid, g); CHKERRQ(nc_check(stat));
        break;
      default:
        SETERRQ1(1, "NC_VAR `%s' not of type NC_INT or NC_FLOAT.\n", name);
    }

    ierr = VecGetArray2d(vindzero, grid.p->Mx, grid.p->My, 0, 0, &ind); CHKERRQ(ierr);
    
    // netCDF concepts of $\Delta x$ and $\Delta y$
    // We have rescaled the grid early in bootstrapFromFile_netCDF() to match the 
    // physical extent of the netCDF file.
    const float ncdx = 2 * grid.p->Lx / (M - 1);
    const float ncdy = 2 * grid.p->Ly / (N - 1);

    for (PetscInt i=0; i < grid.p->Mx; i++) {
      for (PetscInt j=0; j < grid.p->My; j++) {
        const float x = grid.p->dx * (i - grid.p->Mx/2);
        const float y = grid.p->dy * (j - grid.p->My/2);
        if (PetscAbs(x) > grid.p->Lx) {
          SETERRQ1(2, "ncVarToDaVec: x=%f not in bounds.  Grid corrupted.\n", x);
        }
        if (PetscAbs(y) > grid.p->Ly) {
          SETERRQ1(3, "ncVarToDaVec: y=%f not in bounds.  Grid corrupted.\n", y);
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
        if (g != NULL) { // an integer array
          val = w11 * g[i1*N + j1] + w21 * g[i1*N + j2]
                 + w12 * g[i2*N + j1] + w22 * g[i2*N + j2];
          if ((masktool.number_allowed == 1) || (val <= (float)masktool.allowed_levels[0])) {
            val = (float)masktool.allowed_levels[0];
          } else {
            int k=1;        
            while (k < masktool.number_allowed) {
              float mid = ( (float)masktool.allowed_levels[k-1] 
                            + (float)masktool.allowed_levels[k] ) / 2.0;
              if (val < mid) {
                val = (float)masktool.allowed_levels[k-1];
                break;
              }
              k++;
            }
            if (k >= masktool.number_allowed) {
              val = (float)masktool.allowed_levels[k-1];
            }
          }
        } else if (f != NULL) { // a float array
          val = w11 * f[i1*N + j1] + w21 * f[i1*N + j2]
                 + w12 * f[i2*N + j1] + w22 * f[i2*N + j2];
// this way in rev 151 and earlier:
//          val = w11 * f[i1*N + j1] + w12 * f[i1*N + j2]
//                 + w21 * f[i2*N + j1] + w22 * f[i2*N + j2];
        } else {
          SETERRQ(4, "ncvarToDAVec: this should not happen");
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

// bring in the result of applying ncgen.rb to pism_state.cdl (see directory pism/src/netcdf/)
#include "../netcdf/write_attributes.c"

  int s[] = {0, grid.xs, grid.ys, 0};            // Start local block: t dependent
  int c[] = {1, grid.xm, grid.ym, grid.p->Mz};   // Count local block: t dependent
  int cb[] = {1, grid.xm, grid.ym, grid.p->Mbz}; // Count local block: bed

  // Allocate some memory.  We will assume that vectors based on grid.da3 are the largest.
  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
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

  // Allocate some memory.  We will assume that vectors based on grid.da3 are the largest.
  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.p->Mz;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(float), &a_mpi); CHKERRQ(ierr);

// complete the output file as in usual dump
// ("u_id", "v_id", "w_id" defined in included file)
#include "../netcdf/complete_dump.cc"

  // now write additional 3-D diagnostic quantities
  ierr = put_local_var(&grid, ncid, u_id, NC_FLOAT, grid.da3, vu, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, v_id, NC_FLOAT, grid.da3, vv, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = put_local_var(&grid, ncid, w_id, NC_FLOAT, grid.da3, vw, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  
  // We are done with these buffers
  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  if (grid.rank == 0) {
    stat = nc_close(ncid);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  return 0;
}


//! Read a complete saved PISM model state for initialization of an evolution or diagnostic run.
/*! 
When initializing from a NetCDF input file, the file determines 
the number of grid points (Mx,My,Mz,Mbz) and the dimensions of the computational box.   
The user is warned when their command line options "-Mx", "-My", "-Mz", "-Mbz" are overridden.  
 */
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
  // 3-D model quantities
  ierr = get_local_var(&grid, ncid, "T", NC_FLOAT, grid.da3, vT, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "Tb", NC_FLOAT, grid.da3b, vTb, g3b,
                       s, cb, 4, a_mpi, max_a_len); CHKERRQ(ierr);
  ierr = get_local_var(&grid, ncid, "age", NC_FLOAT, grid.da3, vtau, g3,
                       s, c, 4, a_mpi, max_a_len); CHKERRQ(ierr);

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

  ierr = reconfigure_legacy_Mbz(); CHKERRQ(ierr);

  initialized_p = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceModel::readShelfStreamBCFromFile_netCDF(const char *fname) {
  PetscErrorCode  ierr;
  Vec             vbcflag;

  // determine if variables exist in file
  int ncid, stat, 
      maskid, ubarid, vbarid, bcflagid, 
      maskExists=0, ubarExists=0, vbarExists=0, bcflagExists=0;
  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); CHKERRQ(nc_check(stat));
    stat = nc_inq_varid(ncid, "mask", &maskid); maskExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "ubar", &ubarid); ubarExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "vbar", &vbarid); vbarExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "bcflag", &bcflagid); bcflagExists = stat == NC_NOERR;
  }
  ierr = MPI_Bcast(&maskExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ubarExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&vbarExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&bcflagExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  
  if ((ubarExists != 1) || (vbarExists != 1)) {
    SETERRQ1(1,"-shelfstreamBC set but (ubar,vbar) not found in file %s\n",fname);
  }
  if (bcflagExists != 1) {
    SETERRQ1(1,
    "-shelfstreamBC set but bcflag (location of Dirichlet b.c.) not found in file %s\n",
    fname);
  }
  ierr = VecDuplicate(vh, &vbcflag); CHKERRQ(ierr);

  // put read ubar,vbar into da2 Vecs
  // re VecScatter: compare comments in bootstrapFromFile_netCDF()
  Vec vzero;
  VecScatter ctx;
  ierr = VecScatterCreateToZero(g2, &ctx, &vzero); CHKERRQ(ierr);  
  ierr = getIndZero(grid.da2, g2, vzero, ctx); CHKERRQ(ierr);

  if (maskExists) {
    MaskInterp masklevs;
    masklevs.number_allowed = 2;
    masklevs.allowed_levels[0] = MASK_SHEET;
    masklevs.allowed_levels[1] = MASK_FLOATING;
    ierr = ncVarToDAVec(ncid, maskid, grid.da2, vMask, g2, vzero, masklevs); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  mask not found.  Using current values.\n");
               CHKERRQ(ierr);
  }
  ierr = ncVarToDAVec(ncid, ubarid, grid.da2, vubar, g2, vzero); CHKERRQ(ierr);
  ierr = ncVarToDAVec(ncid, vbarid, grid.da2, vvbar, g2, vzero); CHKERRQ(ierr);
  MaskInterp bclevs;
  bclevs.number_allowed = 2;
  bclevs.allowed_levels[0] = 0;
  bclevs.allowed_levels[1] = 1;
  ierr = ncVarToDAVec(ncid, bcflagid, grid.da2, vbcflag, g2, vzero, bclevs); CHKERRQ(ierr);

  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);
  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }

  // now use values in vubar,vvbar, not equal to missing_value, to set boundary conditions by
  // setting corresponding locations to MASK_SHEET and setting vuvbar[0|1] appropriately;
  // set boundary condition which will apply to finite difference system:
  //    staggered grid velocities at MASK_SHEET points which neighbor MASK_FLOATING points
  ierr = VecSet(vuvbar[0],0.0); CHKERRQ(ierr);
  ierr = VecSet(vuvbar[1],0.0); CHKERRQ(ierr);
  PetscScalar **mask, **bc, **ubar, **vbar, **uvbar[2];
  ierr = DAVecGetArray(grid.da2, vbcflag, &bc); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (PetscAbs(bc[i][j] - 1.0) < 0.1) {
        // assume it is really a boundary condition location
        uvbar[0][i-1][j] = ubar[i][j];
        uvbar[0][i][j] = ubar[i][j];
        uvbar[1][i][j-1] = vbar[i][j];
        uvbar[1][i][j] = vbar[i][j];
        mask[i][j] = MASK_SHEET;  // assure that shelf/stream equations not active at this point
      } else {
        uvbar[0][i-1][j] = 0.0;
        uvbar[0][i][j] = 0.0;
        uvbar[1][i][j-1] = 0.0;
        uvbar[1][i][j] = 0.0;        
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vbcflag, &bc); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);    
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);   

  ierr = VecDestroy(vbcflag); CHKERRQ(ierr);

  // update viewers
  ierr = updateViewers(); CHKERRQ(ierr);
  
  // reset initial velocities in shelf for iteration
  ierr = VecSet(vubar,0.0); CHKERRQ(ierr);
  ierr = VecSet(vvbar,0.0); CHKERRQ(ierr);
  return 0;
}

