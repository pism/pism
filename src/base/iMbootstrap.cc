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


//! Read file and use heuristics to initialize PISM from typical 2d data available through remote sensing.
/*! 
This procedure is called when option <tt>-bif</tt> is used.

See chapter 4 of the User's Manual.  We read only 2D information from the bootstrap file.
 */
PetscErrorCode IceModel::bootstrapFromFile_netCDF(const char *fname) {
  PetscErrorCode  ierr;

  const PetscScalar DEFAULT_H_VALUE_NO_VAR = 0.0,  // m
                    DEFAULT_BED_VALUE_NO_VAR = 1.0, // m;  grounded if no bed topo
                    DEFAULT_ACCUM_VALUE_NO_VAR = -0.5 / secpera, // m/s
                    DEFAULT_SURF_TEMP_VALUE_NO_VAR = 263.15, // K
                    DEFAULT_GEOTHERMAL_FLUX_VALUE_NO_VAR = 0.042, // J/m^2 s
                    DEFAULT_UPLIFT_VALUE_NO_VAR = 0.0, // m/s
                    DEFAULT_HMELT_VALUE_NO_VAR = 0.0, // m
                    DEFAULT_TILL_PHI_VALUE_NO_VAR = 15.0; // degrees; tends not to slip
                                                          //   (if -ssa -plastic, which is
                                                          //    not default anyway)

  // start by checking that we have a NetCDF file
  int ncid,  stat;
  PetscInt fileExists=0;
  if (fname == NULL) {
    SETERRQ(1, "no file name given for bootstrapping\n");
  }
  if (grid.rank == 0) {
    stat = nc_open(fname, 0, &ncid); fileExists = (stat == NC_NOERR);
  }
  ierr = MPI_Bcast(&fileExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);  
  if (fileExists) {
    ierr = verbPrintf(2, grid.com, 
       "bootstrapping by PISM default method from file %s\n",fname); 
       CHKERRQ(ierr);
  } else {
    SETERRQ1(2,"bootstrapping file '%s' does not exist\n",fname);
  }
  
  // now allocate our model; when bootstrapFromFile_netCDF() is called 
  //   the grid dimensions (i.e. Mx,My,Mz,Mbz) must already be set
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);

  // determine if dimensions and variables exist in bootstrapping file
  PetscInt tdimExists = 0,
      xdimExists = 0, ydimExists = 0, zdimExists = 0, zbdimExists = 0,
      psExists=0, lonExists=0, latExists=0, accumExists=0, hExists=0, 
      HExists=0, bExists=0, maskExists=0, tillphiExists=0,
      TsExists=0, ghfExists=0, upliftExists=0,
      xExists=0, yExists=0, HmeltExists=0, tExists=0;
  int v_ps, v_lon, v_lat, v_accum, v_h, v_H, v_bed, v_Ts, v_ghf, v_uplift,
      v_x, v_y, v_Hmelt, v_t, v_mask, v_tillphi,
      tdimid, xdimid, ydimid, zdimid, zbdimid;
  if (grid.rank == 0) {
    // use nc_inq_varid to determine whether variable exists in bootstrap file
    // in most cases, the varid itself is discarded
    stat = nc_inq_dimid(ncid, "t", &tdimid); tdimExists = (stat == NC_NOERR);
    stat = nc_inq_dimid(ncid, "x", &xdimid); xdimExists = (stat == NC_NOERR);
    stat = nc_inq_dimid(ncid, "y", &ydimid); ydimExists = (stat == NC_NOERR);
    stat = nc_inq_dimid(ncid, "z", &zdimid); zdimExists = (stat == NC_NOERR);
    stat = nc_inq_dimid(ncid, "zb", &zbdimid); zbdimExists = (stat == NC_NOERR);
    stat = nc_inq_varid(ncid, "x", &v_x); xExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "y", &v_y); yExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "t", &v_t); tExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "polar_stereographic", &v_ps); psExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "lon", &v_lon); lonExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "lat", &v_lat); latExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "acab", &v_accum); accumExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "usurf", &v_h); hExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "thk", &v_H); HExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "topg", &v_bed); bExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "mask", &v_mask); maskExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "artm", &v_Ts); TsExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "bheatflx", &v_ghf); ghfExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "tillphi", &v_tillphi); tillphiExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "dbdt", &v_uplift); upliftExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "bwat", &v_Hmelt); HmeltExists = stat == NC_NOERR;
  }
  // broadcast the existence flags
  ierr = MPI_Bcast(&tdimExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&xdimExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ydimExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&zdimExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&zbdimExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
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
  ierr = MPI_Bcast(&maskExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&TsExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ghfExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&tillphiExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&upliftExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&HmeltExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);

  // if the horizontal dimensions and time (FIXME) are absent then we can not proceed
  if (!tdimExists) {
    SETERRQ1(3,"bootstrapping file '%s' has no time dimension 't'\n",fname);
  }
  if (!xdimExists) {
    SETERRQ1(4,"bootstrapping file '%s' has no horizontal dimension 'x'\n",fname);
  }
  if (!ydimExists) {
    SETERRQ1(5,"bootstrapping file '%s' has no horizontal dimension 'y'\n",fname);
  }
  if (!tExists) {
    SETERRQ1(6,"bootstrapping file '%s' has no variable 't' (=t(t))\n",fname);
  }
  if (!xExists) {
    SETERRQ1(7,"bootstrapping file '%s' has no variable 'x' (=x(x))\n",fname);
  }
  if (!yExists) {
    SETERRQ1(8,"bootstrapping file '%s' has no variable 'y' (=y(y))\n",fname);
  }
 
  // our goal is to create "local interpolation context" from dimensions, 
  //   limits, and lengths
  //   extracted from bootstrap file and from information about the part of the 
  //   grid owned by this processor; this code fills dim[0..4] and bdy[0..6]
  size_t dim[5];  // dimensions in bootstrap NetCDF file
  double bdy[7];  // limits and lengths for bootstrap NetCDF file
  double *z_bif, *zb_bif;
  if ((zdimExists) && (zbdimExists)) {
    ierr = verbPrintf(2, grid.com, 
         "  all dimensions t,x,y,z,zb found in file\n"); CHKERRQ(ierr);
    ierr = nct.get_dims_limits_lengths(ncid, dim, bdy); CHKERRQ(ierr);
    z_bif = new double[dim[3]];
    zb_bif = new double[dim[4]];
    ierr = nct.get_vertical_dims(ncid, dim[3], dim[4], z_bif, zb_bif);
             CHKERRQ(ierr);
  } else if ((!zdimExists) && (!zbdimExists)) {
    ierr = verbPrintf(2, grid.com, 
         "  dimensions t,x,y found in file, but no vertical dimension (z,zb)\n");
         CHKERRQ(ierr);
    ierr = nct.get_dims_limits_lengths_2d(ncid, dim, bdy); CHKERRQ(ierr);
    dim[3] = 1; 
    dim[4] = 1;
    bdy[5] = 0.0;
    bdy[6] = 0.0;
    MPI_Bcast(dim, 5, MPI_LONG, 0, grid.com);
    MPI_Bcast(bdy, 7, MPI_DOUBLE, 0, grid.com);
    z_bif = new double[dim[3]];
    zb_bif = new double[dim[4]];
    z_bif[0] = 0.0;
    zb_bif[0] = 0.0;
  } else {
    SETERRQ(999,"FIXME: need to handle cases of not having ONE of z and zb");
  }
  
  // set year from read-in time variable
  grid.year = bdy[0] / secpera;
  ierr = verbPrintf(2, grid.com, 
          "  time t = %5.4f years found; setting current year\n",
          grid.year); CHKERRQ(ierr);
  ierr = setStartRunEndYearsFromOptions(PETSC_FALSE); CHKERRQ(ierr);

  // runtime options take precedence in setting of -Lx,-Ly,-Lz *including*
  // if initialization is from an input file
  PetscScalar x_scale = (bdy[2] - bdy[1]) / 2.0,
              y_scale = (bdy[4] - bdy[3]) / 2.0,
              z_scale = bdy[6];
  PetscTruth LxSet, LySet, LzSet;
  PetscScalar  x_scale_in, y_scale_in, z_scale_in;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lx", &x_scale_in, &LxSet); CHKERRQ(ierr);
  if (LxSet == PETSC_TRUE)   x_scale = x_scale_in * 1000.0;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Ly", &y_scale_in, &LySet); CHKERRQ(ierr);
  if (LySet == PETSC_TRUE)   y_scale = y_scale_in * 1000.0;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lz", &z_scale_in, &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_TRUE) {
    z_scale = z_scale_in;
  } else {
    ierr = verbPrintf(2, grid.com, 
      "  WARNING: option -Lz should be used to set vertical if bootstrapping ...\n");
      CHKERRQ(ierr);
  }

  // report on resulting computational box, rescale grid, actually create
  //   local interpolation context    
  ierr = verbPrintf(2, grid.com, 
         "  rescaling computational box *for ice* from defaults, -bif file, and\n"
         "    user options to dimensions:\n"
         "    [-%6.2f km, %6.2f km] x [-%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
         x_scale/1000.0,x_scale/1000.0,y_scale/1000.0,y_scale/1000.0,z_scale); 
         CHKERRQ(ierr);
  ierr = determineSpacingTypeFromOptions(PETSC_TRUE); CHKERRQ(ierr);
  ierr = grid.rescale_and_set_zlevels(x_scale, y_scale, z_scale); CHKERRQ(ierr);
  //DEBUG:  ierr = grid.printVertLevels(2); CHKERRQ(ierr);
  LocalInterpCtx lic(ncid, dim, bdy, z_bif, zb_bif, grid);
  //DEBUG:  ierr = lic.printGrid(grid.com); CHKERRQ(ierr);
  delete z_bif;
  delete zb_bif;

  // if polar_stereographic variable exists, then read its attributes
  if (psExists) {
    PetscInt svlfpExists = 0, lopoExists = 0, spExists = 0;
    if (grid.rank == 0) {
      int dummy;
      stat = nc_inq_attid(ncid, v_ps, "straight_vertical_longitude_from_pole",
                              &dummy);
      svlfpExists = (stat == NC_NOERR);
      stat = nc_inq_attid(ncid,v_ps,"latitude_of_projection_origin",&dummy);
      lopoExists = (stat == NC_NOERR);
      stat = nc_inq_attid(ncid,v_ps,"standard_parallel",&dummy);
      spExists = (stat == NC_NOERR);
    }
    ierr = MPI_Bcast(&svlfpExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&lopoExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&spExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
    if (grid.rank == 0) {
       if (svlfpExists) {
         stat = nc_get_att_double(ncid, v_ps, "straight_vertical_longitude_from_pole",
                    &psParams.svlfp); CHKERRQ(nc_check(stat));
       } 
       if (lopoExists) {
         stat = nc_get_att_double(ncid,v_ps,"latitude_of_projection_origin",
                    &psParams.lopo); CHKERRQ(nc_check(stat));
       }
       if (spExists) {
         stat = nc_get_att_double(ncid,v_ps,"standard_parallel",&psParams.sp);
                    CHKERRQ(nc_check(stat));
       }
    }    
    ierr = MPI_Bcast(&psParams.svlfp, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&psParams.lopo, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&psParams.sp, 1, MPI_DOUBLE, 0, grid.com); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
       "  polar stereographic var found; attributes present: svlfp=%d, lopo=%d, sp=%d\n"
       "     values: svlfp = %6.2f, lopo = %6.2f, sp = %6.2f\n",
       svlfpExists, lopoExists, spExists,
       psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr); 
  } else {
    ierr = verbPrintf(2,grid.com,
       "  polar stereo not found, using defaults: svlfp=%6.2f, lopo=%6.2f, sp=%6.2f\n",
       psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr); 
  }

  // now work through all the 2d variables, regridding if present and otherwise setting
  // to default values appropriately
  if (lonExists) {
    ierr = nct.regrid_local_var("lon", 2, lic, grid.da2, vLongitude, g2, false);
             CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vLongitude,"lon","degrees_east",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
      "  WARNING: longitude 'lon' not found; continuing without setting it\n");
      CHKERRQ(ierr);
  }
  if (latExists) {
    ierr = nct.regrid_local_var("lat", 2, lic, grid.da2, vLatitude, g2, false);
             CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vLatitude,"lat","degrees_north",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
      "  WARNING: latitude 'lat' not found; continuing without setting it\n");
      CHKERRQ(ierr);
  }
  if (accumExists) {
    ierr = nct.regrid_local_var("acab", 2, lic, grid.da2, vAccum, g2, false);
             CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vAccum,"acab","m a-1",secpera); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: accumulation rate 'acab' not found; using default constant %7.2f m/a\n",
       DEFAULT_ACCUM_VALUE_NO_VAR * secpera);  CHKERRQ(ierr);
    ierr = VecSet(vAccum, DEFAULT_ACCUM_VALUE_NO_VAR); CHKERRQ(ierr);
  }
  if (HExists) {
    ierr = nct.regrid_local_var("thk", 2, lic, grid.da2, vH, g2, false);
       CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vH,"thk","m",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: thickness 'thk' not found; using default constant %8.2f m\n",
       DEFAULT_H_VALUE_NO_VAR); CHKERRQ(ierr);
    ierr = VecSet(vH, DEFAULT_H_VALUE_NO_VAR); CHKERRQ(ierr); 
  }
  if (bExists) {
    ierr = nct.regrid_local_var("topg", 2, lic, grid.da2, vbed, g2, false);
       CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vbed,"topg","m",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: bedrock elevation 'topg' not found; using default constant %5.4f m\n",
       DEFAULT_BED_VALUE_NO_VAR);  CHKERRQ(ierr);
    ierr = VecSet(vbed, DEFAULT_BED_VALUE_NO_VAR); CHKERRQ(ierr);
  }
  if (maskExists) {
    ierr = verbPrintf(2, grid.com, 
        "  WARNING: 'mask' found; IGNORING IT!\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2, grid.com, 
       "  determining mask by floatation criterion:  grounded ice and ice-free\n"
       "    land marked as 1, floating ice as 3, ice free ocean as 7\n");
     CHKERRQ(ierr);
  if (hExists) {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2, grid.com, 
       "  determining surface elevation by using usurf = topg + thk where grounded\n"
       "    and floatation criterion where floating\n"); CHKERRQ(ierr);
  ierr = setMaskSurfaceElevation_bootstrap(); CHKERRQ(ierr);
  //ierr = VecWAXPY(vh, 1.0, vbed, vH); CHKERRQ(ierr);  // ONLY APPLIED TO GROUNDED
  if (HmeltExists) {
    ierr = nct.regrid_local_var("bwat", 2, lic, grid.da2, vHmelt, g2, false);
            CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vHmelt,"bwat","m",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
        "  WARNING: effective thickness of basal melt water 'bwat' not found;\n"
        "    using default constant %6.3f m\n",
        DEFAULT_HMELT_VALUE_NO_VAR); CHKERRQ(ierr);
    ierr = VecSet(vHmelt,DEFAULT_HMELT_VALUE_NO_VAR); CHKERRQ(ierr);  
  }
  if (TsExists) {
    ierr = nct.regrid_local_var("artm", 2, lic, grid.da2, vTs, g2, false);
       CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vTs,"artm","K",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
       "  WARNING: surface temperature 'artm' not found; using default constant %.2f K\n",
       DEFAULT_SURF_TEMP_VALUE_NO_VAR); CHKERRQ(ierr);
    ierr = VecSet(vTs, DEFAULT_SURF_TEMP_VALUE_NO_VAR); CHKERRQ(ierr);
  }
  if (tillphiExists) {
    ierr = nct.regrid_local_var("tillphi", 2, lic, grid.da2, vtillphi, g2, false);
            CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vtillphi,"tillphi","degrees",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
        "  WARNING: till friction angle 'tillphi' not found; using default constant %6.3f deg\n",
        DEFAULT_TILL_PHI_VALUE_NO_VAR); CHKERRQ(ierr);
    ierr = VecSet(vtillphi,DEFAULT_TILL_PHI_VALUE_NO_VAR); CHKERRQ(ierr);  
  }
  if (ghfExists) {
    ierr = nct.regrid_local_var("bheatflx", 2, lic, grid.da2, vGhf, g2, false);
            CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vGhf,"bheatflx","W m-2",1.0); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: geothermal flux 'bheatflx' not found; using default constant %6.3f W/m^2\n",
       DEFAULT_GEOTHERMAL_FLUX_VALUE_NO_VAR);  CHKERRQ(ierr);
    ierr = VecSet(vGhf, DEFAULT_GEOTHERMAL_FLUX_VALUE_NO_VAR); CHKERRQ(ierr);
  }
  if (upliftExists) {
    ierr = nct.regrid_local_var("dbdt", 2, lic, grid.da2, vuplift, g2, false);
       CHKERRQ(ierr);
    ierr = reportBIFVarFoundMinMax(vuplift,"dbdt","m a-1",secpera); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  WARNING: uplift rate 'dbdt' not found; using default constant %6.3f m/s\n",
       DEFAULT_UPLIFT_VALUE_NO_VAR);  CHKERRQ(ierr);
    ierr = VecSet(vuplift, DEFAULT_UPLIFT_VALUE_NO_VAR); CHKERRQ(ierr);
  }
  
  // fill in temps at depth in reasonable way using surface temps and Ghf
  ierr = verbPrintf(2, grid.com, 
     "  filling in temperatures at depth using quartic guess\n");
     CHKERRQ(ierr);
  ierr = putTempAtDepth(); CHKERRQ(ierr);

  setInitialAgeYears(initial_age_years_default);

  if (grid.rank == 0) {
    stat = nc_close(ncid); CHKERRQ(nc_check(stat));
  }

  ierr = verbPrintf(2, grid.com, "done reading %s; bootstrapping done\n",fname); CHKERRQ(ierr);
  initialized_p = PETSC_TRUE;
  return 0;
}


//! Give a standard report when a variable is found in the -bif file.
PetscErrorCode IceModel::reportBIFVarFoundMinMax(Vec myvar, const char *varname, 
                                 const char *varunits, const PetscScalar factor) {
  PetscErrorCode ierr;
  PetscScalar varmax, varmin, gvarmin, gvarmax;
  
  ierr = VecMin(myvar,PETSC_NULL,&varmin); CHKERRQ(ierr);
  ierr = PetscGlobalMin(&varmin, &gvarmin, grid.com); CHKERRQ(ierr);
  ierr = VecMax(myvar,PETSC_NULL,&varmax); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&varmax, &gvarmax, grid.com); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
    "  variable %s found and regridded; min,max = %8.3f,%8.3f (%s)\n",
    varname,gvarmin*factor,gvarmax*factor,varunits); CHKERRQ(ierr);
  return 0;
}


//! Read certain boundary conditions from a NetCDF file, for diagnostic SSA calculations.
/*!
This is not really a bootstrap procedure, but it has to go somewhere.

For now it is \e only called using "pismd -ross".
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
    SETERRQ1(1,"-ssaBC set but (ubar,vbar) not found in file %s\n",fname);
  }
  if (bcflagExists != 1) {
    SETERRQ1(1,
    "-ssaBC set but bcflag (location of Dirichlet b.c.) not found in file %s\n",
    fname);
  }
  ierr = VecDuplicate(vh, &vbcflag); CHKERRQ(ierr);

  // create "local interpolation context" from dimensions, limits, and lengths extracted from
  //    file and from information about the part of the grid owned by this processor
  size_t dim[5];  // dimensions in NetCDF file
  double bdy[7];  // limits and lengths in NetCDF file
  ierr = nct.get_dims_limits_lengths_2d(ncid, dim, bdy); CHKERRQ(ierr);  // see nc_util.cc
  dim[3] = grid.Mz;  // we ignor any 3D vars, if present, in NetCDF file, and use current grid info
  dim[4] = grid.Mbz;  
  // destructor is called at exit from readShelfStreamBCFromFile_netCDF():
  LocalInterpCtx lic(ncid, dim, bdy, grid.zlevels, grid.zblevels, grid);  

  if (maskExists) {
    MaskInterp masklevs;
    masklevs.number_allowed = 2;
    masklevs.allowed_levels[0] = MASK_SHEET;
    masklevs.allowed_levels[1] = MASK_FLOATING;
    ierr = nct.set_MaskInterp(&masklevs); CHKERRQ(ierr);
    ierr = nct.regrid_local_var("mask", 2, lic, grid.da2, vMask, g2, true); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  mask not found; leaving current values alone ...\n");
               CHKERRQ(ierr);
  }
  ierr = nct.regrid_local_var("ubar", 2, lic, grid.da2, vubar, g2, false); CHKERRQ(ierr);
  ierr = nct.regrid_local_var("vbar", 2, lic, grid.da2, vvbar, g2, false); CHKERRQ(ierr);

  // we have already checked if "bcflag" exists, so just read it
  MaskInterp bclevs;
  bclevs.number_allowed = 2;
  bclevs.allowed_levels[0] = 0;
  bclevs.allowed_levels[1] = 1;
  ierr = nct.set_MaskInterp(&bclevs); CHKERRQ(ierr);
  ierr = nct.regrid_local_var("bcflag", 2, lic, grid.da2, vbcflag, g2, true); CHKERRQ(ierr);

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
        SETERRQ3(1,"Thickness H=%5.4f is negative at point i=%d, j=%d",H[i][j],i,j);
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
        const PetscScalar 
           hgrounded = bed[i][j] + H[i][j],
           hfloating = seaLevel + (1.0 - ice->rho/ocean.rho) * H[i][j];
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
      const PetscScalar beta = (4.0/21.0) * (g / (2.0 * ice->k * HH * HH * HH));
      const PetscScalar alpha = (g / (2.0 * HH * ice->k)) - 2.0 * HH * HH * beta;
      for (PetscInt k = 0; k < ks; k++) {
        const PetscScalar depth = HH - grid.zlevels[k];
        const PetscScalar Tpmp = ice->meltingTemp - ice->beta_CC_grad * depth;
        const PetscScalar d2 = depth * depth;
        T[k] = PetscMin(Tpmp,Ts[i][j] + alpha * d2 + beta * d2 * d2);
      }
      for (PetscInt k = ks; k < grid.Mz; k++) // above ice
        T[k] = Ts[i][j];
      ierr = T3.setValColumnPL(i,j,grid.Mz,grid.zlevels,T); CHKERRQ(ierr);
      
      // set temp within bedrock; if floating then top of bedrock sees ocean,
      //   otherwise it sees the temperature of the base of the ice
      const PetscScalar floating_base = - (ice->rho/ocean.rho) * H[i][j];
      const PetscScalar T_top_bed = (b[i][j] < floating_base)
                                         ? ice->meltingTemp : T[0];
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
This procedure sets the temperatures in the bedrock that would be correct
for our model in steady state.  In steady state there would be a temperature 
at the top of the bed and a flux condition at the bottom
and the temperatures would be linear in between.

Call <tt>Tb3.needAccessToVals()</tt> before and 
<tt>Tb3.doneAccessToVals()</tt> after this routine.
 */
PetscErrorCode IceModel::bootstrapSetBedrockColumnTemp(const PetscInt i, const PetscInt j,
                            const PetscScalar Ttopbedrock, const PetscScalar geothermflux) {
  PetscScalar *Tb;
  Tb = new PetscScalar[grid.Mbz];
  for (PetscInt kb = 0; kb < grid.Mbz; kb++)
    Tb[kb] = Ttopbedrock - (geothermflux / bed_thermal.k) * grid.zblevels[kb];
  PetscErrorCode ierr = Tb3.setInternalColumn(i,j,Tb); CHKERRQ(ierr);
  delete [] Tb;
  return 0;
}

