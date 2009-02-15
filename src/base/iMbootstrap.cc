// Copyright (C) 2004-2009 Jed Brown, Nathan Shemonski, Ed Bueler and
// Constantine Khroulev
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
This procedure is called when option <tt>-boot_from</tt> is used.

See chapter 4 of the User's Manual.  We read only 2D information from the bootstrap file.
 */
PetscErrorCode IceModel::bootstrapFromFile(const char *filename) {
  PetscErrorCode  ierr;

  const PetscScalar DEFAULT_H_VALUE_NO_VAR = 0.0,  // m
                    DEFAULT_BED_VALUE_NO_VAR = 1.0, // m;  grounded if no bed topo
                    DEFAULT_GEOTHERMAL_FLUX_VALUE_NO_VAR = 0.042, // J/m^2 s
                    DEFAULT_UPLIFT_VALUE_NO_VAR = 0.0, // m/s
                    DEFAULT_HMELT_VALUE_NO_VAR = 0.0, // m
                    DEFAULT_TILL_PHI_VALUE_NO_VAR = 15.0; // degrees; tends not to slip
                                                          //   (if -ssa -plastic, which is
                                                          //    not default anyway)

  // start by checking that we have a NetCDF file
  bool file_exists=false;
  NCTool nc(&grid);
  if (filename == NULL) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: No file name given for bootstrapping.\n");
       CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);

  if (file_exists) {
    ierr = verbPrintf(2, grid.com, 
       "bootstrapping by PISM default method from file %s\n",filename); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com,"PISM ERROR: Can't open bootstrapping file '%s'.\n",filename);
       CHKERRQ(ierr);
    PetscEnd();
  }
  
  // now allocate our model; when bootstrapFromFile() is called 
  //   the grid dimensions (i.e. Mx,My,Mz,Mbz) must already be set
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);

  // determine if dimensions and variables exist in bootstrapping file
  bool xdimExists = false, ydimExists = false, zdimExists = false, zbdimExists = false,
       xExists=false, yExists=false, tExists=false,
       hExists=false, maskExists=false;

  ierr = nc.find_dimension("x", NULL, xdimExists); CHKERRQ(ierr);
  ierr = nc.find_dimension("y", NULL, ydimExists); CHKERRQ(ierr);
  ierr = nc.find_dimension("z", NULL, zdimExists); CHKERRQ(ierr);
  ierr = nc.find_dimension("zb", NULL, zbdimExists); CHKERRQ(ierr);
  ierr = nc.find_variable("t", NULL, NULL, tExists); CHKERRQ(ierr);
  ierr = nc.find_variable("x", NULL, NULL, xExists); CHKERRQ(ierr);
  ierr = nc.find_variable("y", NULL, NULL, yExists); CHKERRQ(ierr);
  ierr = nc.find_variable("usurf", "surface_altitude", NULL,  hExists); CHKERRQ(ierr);
  ierr = nc.find_variable("mask", NULL, NULL, maskExists); CHKERRQ(ierr);

  // if the horizontal dimensions are absent then we can not proceed
  if (!xdimExists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'x'\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (!ydimExists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'y'\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (!xExists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no variable 'x' (=x(x))\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (!yExists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no variable 'y' (=y(y))\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
 
  // our goal is to create "local interpolation context" from dimensions,
  //   limits, and lengths extracted from bootstrap file and from information
  //   about the part of the grid owned by this processor
  grid_info g;
  ierr = nc.get_grid_info(g); CHKERRQ(ierr); // g.z_max is set to 0 if z does
					     // not exist

  if (tExists) {
    grid.year = g.time / secpera; // set year from read-in time variable
    ierr = verbPrintf(2, grid.com, 
		      "  time t = %5.4f years found; setting current year\n",
		      grid.year); CHKERRQ(ierr);
  } else {
    grid.year = 0.0;
    ierr = verbPrintf(2, grid.com, 
		      "  time dimension was not found; setting current year to t = 0.0 years\n",
		      grid.year); CHKERRQ(ierr);
  }
  ierr = setStartRunEndYearsFromOptions(PETSC_TRUE); CHKERRQ(ierr);

  // Record the coordinates of the grid center:
  grid.x0 = g.x0;
  grid.y0 = g.y0;

  // runtime options take precedence in setting of -Lx,-Ly,-Lz *including*
  // if initialization is from an input file
  PetscScalar x_scale = (g.x_max - g.x_min) / 2.0,
              y_scale = (g.y_max - g.y_min) / 2.0,
              z_scale = g.z_max;
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
         "  rescaling computational box for ice from defaults, -boot_from file, and\n"
         "    user options to dimensions:\n"
         "    [-%6.2f km, %6.2f km] x [-%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
         x_scale/1000.0,x_scale/1000.0,y_scale/1000.0,y_scale/1000.0,z_scale); 
         CHKERRQ(ierr);
  ierr = determineSpacingTypeFromOptions(PETSC_TRUE); CHKERRQ(ierr);
  ierr = grid.rescale_and_set_zlevels(x_scale, y_scale, z_scale); CHKERRQ(ierr);
  
  // now we have enough to actually create the lic
  // IceModel::bootstrapLIC is now a valid pointer which can be reused to
  // get more info out of the the -boot_from file
  bootstrapLIC = new LocalInterpCtx(g, NULL, NULL, grid);

  ierr = nc.read_polar_stereographic(psParams.svlfp,
				     psParams.lopo,
				     psParams.sp,
				     true); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  // now work through all the 2d variables, regridding if present and otherwise setting
  // to default values appropriately

  if (maskExists) {
    ierr = verbPrintf(2, grid.com, 
		      "  WARNING: 'mask' found; IGNORING IT!\n"); CHKERRQ(ierr);
  }
  if (hExists) {
    ierr = verbPrintf(2, grid.com, 
		      "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n"); CHKERRQ(ierr);
  }


  ierr = verbPrintf(2, grid.com, 
		    "processing 2D model state variables...\n"); CHKERRQ(ierr);

  ierr = vLongitude.regrid(filename, *bootstrapLIC, false); CHKERRQ(ierr);
  ierr =  vLatitude.regrid(filename, *bootstrapLIC, false); CHKERRQ(ierr);
  ierr =         vH.regrid(filename, *bootstrapLIC, DEFAULT_H_VALUE_NO_VAR); CHKERRQ(ierr);
  ierr =       vbed.regrid(filename, *bootstrapLIC, DEFAULT_BED_VALUE_NO_VAR); CHKERRQ(ierr);
  ierr =     vHmelt.regrid(filename, *bootstrapLIC, DEFAULT_HMELT_VALUE_NO_VAR); CHKERRQ(ierr);
  ierr =   vtillphi.regrid(filename, *bootstrapLIC, DEFAULT_TILL_PHI_VALUE_NO_VAR); CHKERRQ(ierr);
  ierr =       vGhf.regrid(filename, *bootstrapLIC, DEFAULT_GEOTHERMAL_FLUX_VALUE_NO_VAR); CHKERRQ(ierr);
  ierr =    vuplift.regrid(filename, *bootstrapLIC, DEFAULT_UPLIFT_VALUE_NO_VAR); CHKERRQ(ierr);

  // set mask and h; tell user what happened:
  ierr = setMaskSurfaceElevation_bootstrap(); CHKERRQ(ierr);

  setInitialAgeYears(initial_age_years_default);
  
  // fill in temps at depth in reasonable way using surface temps (from atmosPCC) and Ghf
  ierr = verbPrintf(2, grid.com, 
     "  filling in temperatures at depth using quartic guess\n");
     CHKERRQ(ierr);
#if 0
#else
  // FIXME: kludge to deal with hosed initialization sequence
  ierr = vWork2d[0].set_name("artm"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs(
            "climate_state", "temperature at ice surface but below firn",
            "K", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].regrid(filename, *bootstrapLIC, 273.15); CHKERRQ(ierr);
#endif
  ierr = putTempAtDepth(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "done reading %s; bootstrapping done\n",filename); CHKERRQ(ierr);
  initialized_p = PETSC_TRUE;
  return 0;
}


//! Read certain boundary conditions from a NetCDF file, for diagnostic SSA calculations.
/*!
This is not really a bootstrap procedure, but it has to go somewhere.

For now it is \e only called using "pismd -ross".
 */
PetscErrorCode IceModel::readShelfStreamBCFromFile(const char *filename) {
  PetscErrorCode  ierr;
  IceModelVec2 vbcflag;
  NCTool nc(&grid);
  bool file_exists;

  // determine if variables exist in file
  int maskid, ubarid, vbarid, bcflagid;

  bool maskExists=false, ubarExists=false, vbarExists=false, bcflagExists=false; 

  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Can't open file '%s'.\n", filename);
    CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = nc.find_variable("mask", NULL, &maskid, maskExists); CHKERRQ(ierr);
  ierr = nc.find_variable("ubar", NULL, &ubarid, ubarExists); CHKERRQ(ierr);
  ierr = nc.find_variable("vbar", NULL, &vbarid, vbarExists); CHKERRQ(ierr);
  ierr = nc.find_variable("bcflag", NULL, &bcflagid, bcflagExists); CHKERRQ(ierr);
  
  if ((ubarExists != 1) || (vbarExists != 1)) {
    ierr = PetscPrintf(grid.com,"-ssaBC set but (ubar,vbar) not found in file %s\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (bcflagExists != 1) {
    ierr = PetscPrintf(grid.com,
		       "-ssaBC set but bcflag (location of Dirichlet b.c.) not found in file %s\n",
		       filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = vbcflag.create(grid, "bcflag", true); CHKERRQ(ierr);

  // create "local interpolation context" from dimensions, limits, and lengths extracted from
  //    file and from information about the part of the grid owned by this processor
  grid_info g;
  ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);  // see nc_util.cc
  // destructor is called at exit from readShelfStreamBCFromFile():
  LocalInterpCtx lic(g, NULL, NULL, grid); // 2D only

  if (maskExists) {
    vMask.interpolation_mask.number_allowed = 2;
    vMask.interpolation_mask.allowed_levels[0] = MASK_SHEET;
    vMask.interpolation_mask.allowed_levels[1] = MASK_FLOATING;
    vMask.use_interpolation_mask = true;
    ierr = vMask.regrid(filename, lic, true); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(3, grid.com, "  mask not found; leaving current values alone ...\n");
               CHKERRQ(ierr);
  }
  ierr = vubar.regrid(filename, lic, true); CHKERRQ(ierr);
  ierr = vvbar.regrid(filename, lic, true); CHKERRQ(ierr);

  // we have already checked if "bcflag" exists, so just read it
  vbcflag.interpolation_mask.number_allowed = 2;
  vbcflag.interpolation_mask.allowed_levels[0] = 0;
  vbcflag.interpolation_mask.allowed_levels[1] = 1;
  vbcflag.use_interpolation_mask = true;
  ierr = vbcflag.regrid(filename, lic, true); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  // now use values in vubar,vvbar, not equal to missing_value, to set boundary conditions by
  // setting corresponding locations to MASK_SHEET and setting vuvbar[0|1] appropriately;
  // set boundary condition which will apply to finite difference system:
  //    staggered grid velocities at MASK_SHEET points which neighbor MASK_FLOATING points
  ierr = vuvbar[0].set(0.0); CHKERRQ(ierr);
  ierr = vuvbar[1].set(0.0); CHKERRQ(ierr);
  PetscScalar **mask, **bc, **ubar, **vbar, **uvbar[2];
  ierr = vbcflag.get_array(bc); CHKERRQ(ierr);    
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vubar.get_array(ubar); CHKERRQ(ierr);
  ierr = vvbar.get_array(vbar); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);

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
  ierr =   vbcflag.end_access(); CHKERRQ(ierr);    
  ierr =     vMask.end_access(); CHKERRQ(ierr);
  ierr =     vubar.end_access(); CHKERRQ(ierr);
  ierr =     vvbar.end_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);

  ierr = vbcflag.destroy(); CHKERRQ(ierr);

  // update viewers
  ierr = updateViewers(); CHKERRQ(ierr);
  return 0;
}


//! Determine surface and mask according to information in bootstrap file and options.
PetscErrorCode IceModel::setMaskSurfaceElevation_bootstrap() {
    PetscErrorCode ierr;
  PetscScalar **h, **bed, **H, **mask;

  ierr = verbPrintf(2, grid.com, 
    "  determining surface elevation by  usurf = topg + thk  where grounded\n"
    "    and by floatation crit  usurf = (1-rho_i/rho_w) thk  where floating\n"); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
           "  determining mask:  grounded ice and ice-free land marked as 1,\n"); CHKERRQ(ierr);
  if (doOceanKill == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com,
           "    floating ice as 3, and ice free ocean as 7 (from -ocean_kill)\n"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
           "    floating ice and ice-free ocean as 3\n"); CHKERRQ(ierr);
  }

  if (oceanPCC == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");  }
  PetscReal currentSeaLevel = oceanPCC->reportSeaLevelElevation();
           
  ierr = vh.get_array(h); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vbed.get_array(bed); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // take this opportunity to check that H[i][j] >= 0
      if (H[i][j] < 0.0) {
        SETERRQ3(2,"Thickness H=%5.4f is negative at point i=%d, j=%d",H[i][j],i,j);
      }
      
      if (H[i][j] < 0.001) {  // if no ice
        if (bed[i][j] < 0.0) {
          h[i][j] = 0.0;
          mask[i][j] = (doOceanKill == PETSC_TRUE) ? MASK_FLOATING_OCEAN0 : MASK_FLOATING;
        } else {
          h[i][j] = bed[i][j];
          mask[i][j] = MASK_SHEET;
        } 
      } else { // if positive ice thickness then check floatation criterion
        const PetscScalar 
           hgrounded = bed[i][j] + H[i][j],
           hfloating = currentSeaLevel + (1.0 - ice->rho/ocean.rho) * H[i][j];
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

  ierr =    vh.end_access(); CHKERRQ(ierr);
  ierr =    vH.end_access(); CHKERRQ(ierr);
  ierr =  vbed.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  // go ahead and communicate mask and surface elev now; may be redundant communication?
  ierr =    vh.beginGhostComm(); CHKERRQ(ierr);
  ierr =    vh.endGhostComm(); CHKERRQ(ierr);
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
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

  IceModelVec2    *pccTs;
#if 0
  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, dt * secpera, (void*)(&info_atmoscoupler), pccTs); CHKERRQ(ierr);
  } else {  SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");  }
#else
  // FIXME: kludge to deal with hosed initialization sequence
  pccTs = &(vWork2d[0]);
#endif

  ierr = pccTs->get_array(Ts);  CHKERRQ(ierr);
  ierr =   vH.get_array(H);   CHKERRQ(ierr);
  ierr = vbed.get_array(b);   CHKERRQ(ierr);
  ierr = vGhf.get_array(Ghf); CHKERRQ(ierr);

  ierr =  T3.begin_access(); CHKERRQ(ierr);
  ierr = Tb3.begin_access(); CHKERRQ(ierr);
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
  ierr =   vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr =   T3.end_access(); CHKERRQ(ierr);
  ierr =  Tb3.end_access(); CHKERRQ(ierr);
  ierr = pccTs->end_access(); CHKERRQ(ierr);

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

Call <tt>Tb3.begin_access()</tt> before and 
<tt>Tb3.end_access()</tt> after this routine.
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

