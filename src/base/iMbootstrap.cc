// Copyright (C) 2004-2008 Jed Brown and Nathan Shemonski and Ed Bueler
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


//! Create a temperature field within ice and bedrock from given surface temperature and geothermal flux maps.
/*!
In bootstrapping we need to guess about the temperature within the ice and bedrock if surface temperature
and geothermal flux maps are given.  This rule is heuristic but seems to work well anyway.  Full 
bootstrapping will start from the temperature computed by this procedure and then run for a long time 
(e.g. \f$10^5\f$ years), with fixed geometry, to get closer to thermomechanically coupled equilibrium.
See the part of the <i>User's Manual</i> on EISMINT-Greenland.

Consider a horizontal grid point <tt>i,j</tt>.  Suppose the surface temperature \f$T_s\f$ and the geothermal
flux \f$g\f$ are given at that grid point.  Within the corresponding column, denote the temperature
by \f$T(z)\f$ for some elevation \f$z\f$ above the base of the ice.  (Note ice corresponds to \f$z>0\f$ while
bedrock has \f$z<0\f$.)  Apply the rule that \f$T(z)=T_s\f$ is \f$z\f$ is above the top of the ice (at 
\f$z=H\f$).  

Within the ice, set
	\f[T(z) = T_s + \alpha (H-z)^2 + \beta (H-z)^4\f]
where \f$\alpha,\beta\f$ are chosen so that
	\f[\frac{\partial T}{\partial z}\Big|_{z=0} = - \frac{g}{k_i}\f]
and 
   \f[\frac{\partial T}{\partial z}\Big|_{z=H/4} = - \frac{g}{2 k_i}.\f]
The point of the second condition is our observation that, in observed ice, the rate of decrease 
in ice temperature with elevation is significantly decreased at only one quarter of the ice thickness above 
the base.  

The temperature within the ice is not allowed to exceed the pressure-melting temperature.

Note that the above heuristic rule for ice determines \f$T(0)\f$.  Within the bedrock our rule is that 
the rate of change with depth is exactly the geothermal flux:
   \f[T(z) = T(0) - \frac{g}{k_r} z.\f]
Note that \f$z\f$ here is negative, so the temperature increases as one goes down into the bed.
 */
PetscErrorCode IceModel::putTempAtDepth() {
  PetscErrorCode  ierr;
  PetscScalar     **H, **b, **Ts, **Ghf;

  PetscScalar *T;
  T = new PetscScalar[grid.Mz];

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = T3.needAccessToVals(); CHKERRQ(ierr);
  ierr = Tb3.needAccessToVals(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar HH = H[i][j];
      const PetscInt    ks = grid.kBelowHeight(HH);
      
      // within ice
      const PetscScalar g = Ghf[i][j];
      const PetscScalar beta = (4.0/21.0) * (g / (2.0 * ice.k * HH * HH * HH));
      const PetscScalar alpha = (g / (2.0 * HH * ice.k)) - 2.0 * HH * HH * beta;
      for (PetscInt k = 0; k < ks; k++) {
        const PetscScalar depth = HH - grid.zlevels[k];
        const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * depth;
        const PetscScalar d2 = depth * depth;
        T[k] = PetscMin(Tpmp,Ts[i][j] + alpha * d2 + beta * d2 * d2);
      }
      for (PetscInt k = ks; k < grid.Mz; k++) // above ice
        T[k] = Ts[i][j];
      ierr = T3.setValColumnPL(i,j,grid.Mz,grid.zlevels,T); CHKERRQ(ierr);
      
      // set temp within bedrock; if floating then top of bedrock sees ocean,
      //   otherwise it sees the temperature of the base of the ice
      const PetscScalar floating_base = - (ice.rho/ocean.rho) * H[i][j];
      const PetscScalar T_top_bed = (b[i][j] < floating_base)
                                         ? ice.meltingTemp : T[0];
      ierr = bootstrapSetBedrockColumnTemp(i,j,T_top_bed,Ghf[i][j]); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = T3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = Tb3.doneAccessToVals(); CHKERRQ(ierr);

  delete [] T;
  
  ierr = T3.beginGhostComm(); CHKERRQ(ierr);
  ierr = T3.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


//! Set the temperatures in a column of bedrock based on a temperature at the top and a geothermal flux.
/*! 
This procedure sets the temperatures in the bedrock that would be correct for our model in steady state.
In steady state there would be a temperature at the top of the bed and a flux condition at the bottom
and the temperatures would be linear in between.

Call <tt>Tb3.needAccessToVals()</tt> before and <tt>Tb3.doneAccessToVals()</tt> after this routine.
 */
PetscErrorCode IceModel::bootstrapSetBedrockColumnTemp(const PetscInt i, const PetscInt j,
                            const PetscScalar Ttopbedrock, const PetscScalar geothermflux) {
  PetscScalar *Tb;
  Tb = new PetscScalar[grid.Mbz];
  for (PetscInt kb = 0; kb < grid.Mbz; kb++)
    Tb[kb] = Ttopbedrock - (geothermflux / bedrock.k) * grid.zblevels[kb];
  PetscErrorCode ierr = Tb3.setInternalColumn(i,j,Tb); CHKERRQ(ierr);
  delete [] Tb;
  return 0;
}


//! Use heuristics to set up ice sheet fields needed for actual initialization.
/*! 
This procedure is called when option <tt>-bif</tt> is used.

Note this procedure uses var_to_da_vec() which is DEPRECATED.
 */
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
    stat = nc_inq_varid(ncid, "acab", &v_accum); accumExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "usurf", &v_h); hExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "thk", &v_H); HExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "topg", &v_bed); bExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "artm", &v_Ts); TsExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "bheatflx", &v_ghf); ghfExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "dbdt", &v_uplift); upliftExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "balvel", &v_balvel); balvelExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "bwat", &v_Hmelt); HmeltExists = stat == NC_NOERR;
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
    ierr = verbPrintf(2,grid.com,
               "  polar stereographic found: svlfp = %6.2f, lopo = %6.2f, sp = %6.2f\n",
               psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr); 
  } else {
    ierr = verbPrintf(2,grid.com,
               "  polar stereo not found, using defaults: svlfp=%6.2f, lopo=%6.2f, sp=%6.2f\n",
               psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr); 
  }

  // set time (i.e. grid.year, IceModel::startYear, and IceModel::endYear)
  if (tExists) {
    const size_t idx = 0;
    double begin_t;
    if (grid.rank == 0) {
      nc_get_var1_double(ncid, v_t, &idx, &begin_t);
    }
    ierr = MPI_Bcast(&begin_t, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    grid.year = begin_t / secpera;
    ierr = verbPrintf(2, grid.com, 
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
    ierr = nct.get_ends_1d_var(ncid, v_x, &first, &last, grid.com); CHKERRQ(ierr);
    x_scale = (last - first) / 2.0;
  } else {
    ierr = verbPrintf(2, grid.com,"  WARNING: variable x(x) not found in bootstrap file.\n");
             CHKERRQ(ierr);
  }
  if (yExists) {
    ierr = nct.get_ends_1d_var(ncid, v_y, &first, &last, grid.com); CHKERRQ(ierr);
    y_scale = (last - first) / 2.0;
  } else {
    ierr = verbPrintf(2, grid.com,"  WARNING: variable y(y) not found in bootstrap file.\n");
             CHKERRQ(ierr);
  }

  // runtime options take precedence in setting of -Lx,-Ly,-Lz *including*
  // if initialization is from an input file
  PetscTruth LxSet, LySet, LzSet;
  PetscScalar  x_scale_in, y_scale_in, z_scale_in;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lx", &x_scale_in, &LxSet); CHKERRQ(ierr);
  if (LxSet == PETSC_TRUE)   x_scale = x_scale_in * 1000.0;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Ly", &y_scale_in, &LySet); CHKERRQ(ierr);
  if (LySet == PETSC_TRUE)   y_scale = y_scale_in * 1000.0;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lz", &z_scale_in, &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_TRUE)   z_scale = z_scale_in;
    
  ierr = verbPrintf(2, grid.com, 
         "  rescaling computational box *for ice* from defaults, -bif file, and\n"
         "    user options to dimensions:\n"
         "      [-%6.2f km, %6.2f km] x [-%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
         x_scale/1000.0,x_scale/1000.0,y_scale/1000.0,y_scale/1000.0,z_scale); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
         "  resetting vertical levels using user options and grid.rescale_and_set_zlevels()\n");
         CHKERRQ(ierr);
  ierr = determineSpacingTypeFromOptions(); CHKERRQ(ierr);
  ierr = grid.rescale_and_set_zlevels(x_scale, y_scale, z_scale); CHKERRQ(ierr);

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
    ierr = nct.var_to_da_vec(grid, ncid, v_lon, grid.da2, vLongitude, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "  WARNING: longitude lon not found; continuing without it\n");  CHKERRQ(ierr);
  }
  if (latExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_lat, grid.da2, vLatitude, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "  WARNING: latitude lat not found; continuing without it\n");  CHKERRQ(ierr);
  }
  if (accumExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_accum, grid.da2, vAccum, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
               "  WARNING: accumulation rate acab not found; using default %7.2f m/a\n",
               DEFAULT_ACCUM_VALUE_MISSING * secpera);  CHKERRQ(ierr);
    ierr = VecSet(vAccum, DEFAULT_ACCUM_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (HExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_H, grid.da2, vH, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
               "  WARNING: thickness thk not found; using default %8.2f\n",
               DEFAULT_H_VALUE_MISSING); CHKERRQ(ierr);
    ierr = VecSet(vH, DEFAULT_H_VALUE_MISSING); CHKERRQ(ierr); 
  }
  if (bExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_bed, grid.da2, vbed, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
             "  WARNING: bedrock elevation topg not found; using default %f\n",
             DEFAULT_BED_VALUE_MISSING);  CHKERRQ(ierr);
    ierr = VecSet(vbed, DEFAULT_BED_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (hExists) {
    ierr = verbPrintf(2, grid.com, 
          "  WARNING: ignoring values found for surface elevation; using usurf = topg + thk\n");
          CHKERRQ(ierr);
  }
  if (TsExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_Ts, grid.da2, vTs, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
             "  WARNING: surface temperature artm not found; using default %7.2f K\n",
             DEFAULT_SURF_TEMP_VALUE_MISSING); CHKERRQ(ierr);
    ierr = VecSet(vTs, DEFAULT_SURF_TEMP_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (ghfExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_ghf, grid.da2, vGhf, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
             "  WARNING: geothermal flux bheatflx not found; using default %6.3f W/m^2\n",
             DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING);  CHKERRQ(ierr);
    ierr = VecSet(vGhf, DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING); CHKERRQ(ierr);
  }
  if (upliftExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_uplift, grid.da2, vuplift, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: uplift rate dbdt not found; filling with zero\n");  CHKERRQ(ierr);
    ierr = VecSet(vuplift, 0.0); CHKERRQ(ierr);
  }
  if (HmeltExists) {
    ierr = nct.var_to_da_vec(grid, ncid, v_Hmelt, grid.da2, vHmelt, g2, vzero); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
        "  WARNING: effective thickness of basal melt water bwat not found; filling with zero\n"); 
        CHKERRQ(ierr);
    ierr = VecSet(vHmelt,0.0); CHKERRQ(ierr);  
  }
  // done with scatter context
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);
  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }
  ierr = verbPrintf(3, grid.com, "  done reading .nc file\n"); CHKERRQ(ierr);

//  ierr = cleanJustNan_legacy(); CHKERRQ(ierr);  // needed for ant_init.nc; should be harmless for other

  ierr = verbPrintf(2, grid.com, 
            "  determining mask by floating criterion; grounded ice marked as SIA (=1)\n");
            CHKERRQ(ierr);
  ierr = setMaskSurfaceElevation_bootstrap(); CHKERRQ(ierr);
  
  // fill in temps at depth in reasonable way using surface temps and Ghf
  ierr = verbPrintf(2, grid.com, 
             "  filling in temperatures at depth using surface temperatures and quartic guess\n"); CHKERRQ(ierr);
  ierr = putTempAtDepth(); CHKERRQ(ierr);
  // fill in other 3D fields
  setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
  setConstantGrainSize(DEFAULT_GRAIN_SIZE);

  verbPrintf(2, grid.com, "bootstrapping by PISM default method done\n");
  initialized_p = PETSC_TRUE;
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


//! Read certain boundary conditions from a NetCDF file, especially for diagnostic SSA calculations.
/*!
This is not really a bootstrap procedure, but it has to go somewhere.

Note this procedure uses var_to_da_vec() which is DEPRECATED.
 */
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
    ierr = nct.var_to_da_vec(grid, ncid, maskid, grid.da2, vMask, g2, vzero, masklevs); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  mask not found.  Using current values.\n");
               CHKERRQ(ierr);
  }
  ierr = nct.var_to_da_vec(grid, ncid, ubarid, grid.da2, vubar, g2, vzero); CHKERRQ(ierr);
  ierr = nct.var_to_da_vec(grid, ncid, vbarid, grid.da2, vvbar, g2, vzero); CHKERRQ(ierr);
  MaskInterp bclevs;
  bclevs.number_allowed = 2;
  bclevs.allowed_levels[0] = 0;
  bclevs.allowed_levels[1] = 1;
  ierr = nct.var_to_da_vec(grid, ncid, bcflagid, grid.da2, vbcflag, g2, vzero, bclevs); CHKERRQ(ierr);

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

