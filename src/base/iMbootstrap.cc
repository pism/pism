// Copyright (C) 2004-2011 Jed Brown, Nathan Shemonski, Ed Bueler and
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

#include "iceModel.hh"
#include "PISMIO.hh"

//! Read file and use heuristics to initialize PISM from typical 2d data available through remote sensing.
/*! 
This procedure is called by the base class when option <tt>-boot_file</tt> is used.

See chapter 4 of the User's Manual.  We read only 2D information from the bootstrap file.
 */
PetscErrorCode IceModel::bootstrapFromFile(const char *filename) {
  PetscErrorCode  ierr;

  // Bootstrap 2D fields:
  ierr = bootstrap_2d(filename); CHKERRQ(ierr);

  // Regrid 2D fields:
  ierr = regrid(2); CHKERRQ(ierr);

  // Check the consistency of geometry fields:
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); 

  // If ocean_kill is set, mark ice-free ocean with a special mask value:
  ierr = mark_ocean_at_time_0(); CHKERRQ(ierr);

  // Update couplers (because heuristics in bootstrap_3d() might need boundary
  // conditions provided by couplers):
  if (surface != NULL) {
    ierr = surface->update(grid.year, 0); CHKERRQ(ierr);
  } else SETERRQ(1, "surface == NULL");
  if (ocean != NULL) {
    ierr = ocean->update(grid.year, 0); CHKERRQ(ierr);
  } else SETERRQ(1, "ocean == NULL");

  // Fill 3D fields using heuristics:
  ierr = bootstrap_3d(); CHKERRQ(ierr);

  ierr = regrid(3); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "done reading %s; bootstrapping done\n",filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::bootstrap_2d(const char *filename) {
  PetscErrorCode ierr;

  PISMIO nc(&grid);
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
		    "bootstrapping by PISM default method from file %s\n",filename); CHKERRQ(ierr);

  // report on resulting computational box, rescale grid, actually create a
  // local interpolation context
  ierr = verbPrintf(2, grid.com, 
         "  rescaling computational box for ice from -boot_file file and\n"
         "    user options to dimensions:\n"
         "    [-%6.2f km, %6.2f km] x [-%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
         grid.Lx/1000.0,grid.Lx/1000.0,grid.Ly/1000.0,grid.Ly/1000.0,grid.Lz); 
         CHKERRQ(ierr);

  bool hExists=false, maskExists=false;
  ierr = nc.find_variable("usurf", "surface_altitude", NULL,  hExists); CHKERRQ(ierr);
  ierr = nc.find_variable("mask", NULL, maskExists); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // now work through all the 2d variables, regridding if present and otherwise
  // setting to default values appropriately

  if (maskExists) {
    ierr = verbPrintf(2, grid.com, 
		      "  WARNING: 'mask' found; IGNORING IT!\n"); CHKERRQ(ierr);
  }
  if (hExists) {
    ierr = verbPrintf(2, grid.com, 
		      "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n");
		      CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, 
		    "  reading 2D model state variables by regridding ...\n"); CHKERRQ(ierr);

  ierr = vLongitude.regrid(filename, false); CHKERRQ(ierr);
  ierr =  vLatitude.regrid(filename, false); CHKERRQ(ierr);
  ierr =         vH.regrid(filename, 
                           config.get("bootstrapping_H_value_no_var")); CHKERRQ(ierr);
  ierr =       vbed.regrid(filename,  
                           config.get("bootstrapping_bed_value_no_var")); CHKERRQ(ierr);
  ierr =     vHmelt.regrid(filename,  
                           config.get("bootstrapping_Hmelt_value_no_var")); CHKERRQ(ierr);
  ierr =       vbmr.regrid(filename,  
                           config.get("bootstrapping_bmelt_value_no_var")); CHKERRQ(ierr);
  ierr =   vtillphi.regrid(filename,  
                           config.get("bootstrapping_tillphi_value_no_var")); CHKERRQ(ierr);
  ierr =       vGhf.regrid(filename,  
                           config.get("bootstrapping_geothermal_flux_value_no_var"));
  CHKERRQ(ierr);
  ierr =    vuplift.regrid(filename,  
                           config.get("bootstrapping_uplift_value_no_var")); CHKERRQ(ierr);

  if (config.get_flag("part_grid")) {
    // if part_grid is "on", set fields tracking contents of partially-filled
    // cells to zero. Note that the contents of these fields are
    // grid-dependent, so we don't want to read them from a bootstrapping file
    // using linear interpolation.
    //ierr = vHav.set(0.0); CHKERRQ(ierr);
    ierr = vHref.set(0.0); CHKERRQ(ierr);
    if (config.get_flag("part_redist")) { ierr = vHresidual.set(0.0); CHKERRQ(ierr); }
  }

  if (config.get_flag("kill_icebergs")) {
    // will be updated in updateSurfaceElevationAndMask()
    ierr = vIcebergMask.set(ICEBERGMASK_NOT_SET); CHKERRQ(ierr);
  }

  if (config.get_flag("do_eigen_calving")) {
    // will be updated in updateSurfaceElevationAndMask()
    ierr = vPrinStrain1.set(0.0); CHKERRQ(ierr);
    ierr = vPrinStrain2.set(0.0); CHKERRQ(ierr);
  }

  if (config.get_flag("dirichlet_bc")) {
    ierr = vBCMask.regrid(filename,
                             config.get("bootstrapping_BCMask_value_no_var")); CHKERRQ(ierr);
    ierr = vBCvel.regrid(filename,
                             config.get("bootstrapping_BCvel_value_no_var")); CHKERRQ(ierr);
  }



  bool Lz_set;
  ierr = PISMOptionsIsSet("-Lz", Lz_set); CHKERRQ(ierr);
  if ( !Lz_set ) {
    PetscReal thk_min, thk_max;
    ierr = vH.range(thk_min, thk_max); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
		      "  Setting Lz to 1.5 * max(ice thickness) = %3.3f meters...\n",
		      1.5 * thk_max);


    grid.Lz = 1.5 * thk_max;

    ierr = grid.compute_vertical_levels();

    CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::bootstrap_3d() {
  PetscErrorCode ierr;

  // set the initial age of the ice if appropriate
  if (config.get_flag("do_age")) {
    ierr = verbPrintf(2, grid.com, 
      "  setting initial age to %.4f years\n", config.get("initial_age_of_ice_years"));
      CHKERRQ(ierr);
    tau3.set(config.get("initial_age_of_ice_years") * secpera);
  }
  
  ierr = verbPrintf(2, grid.com, 
     "  filling in ice and bedrock temperatures using surface temps and quartic guess\n");
     CHKERRQ(ierr);
  ierr = putTempAtDepth(); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods") == false) {
    ierr = verbPrintf(2, grid.com,
		      "  ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");
    CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::mark_ocean_at_time_0() {
  PetscErrorCode ierr;

  if (config.get_flag("ocean_kill") == false) return 0;

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vMask(i,j) == MASK_ICE_FREE_OCEAN)
        vMask(i,j) = MASK_OCEAN_AT_TIME_0;
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);

  // Communicate mask values (just in case):
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
bedrock has \f$z<0\f$.)  Apply the rule that \f$T(z)=T_s\f$ is \f$z\f$ is above the top of the ice (at \f$z=H\f$).  

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

FIXME task #7297
 */
PetscErrorCode IceModel::putTempAtDepth() {
  PetscErrorCode  ierr;

  PetscScalar *T = new PetscScalar[grid.Mz];
  const bool do_cold = config.get_flag("do_cold_ice_methods");

  if (surface != NULL) {
    ierr = surface->ice_surface_temperature(artm); CHKERRQ(ierr);
  } else {
    SETERRQ(1, "PISM ERROR: surface == NULL");
  }

  IceModelVec3 *result;
  if (do_cold) 
    result = &T3;
  else
    result = &Enth3;

  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr =   vH.begin_access();   CHKERRQ(ierr);
  ierr = vbed.begin_access();   CHKERRQ(ierr);
  ierr = vGhf.begin_access(); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar HH = vH(i,j);
      const PetscInt    ks = grid.kBelowHeight(HH);
      
      // within ice
      const PetscScalar g = vGhf(i,j);
      const PetscScalar beta = (4.0/21.0) * (g / (2.0 * ice->k * HH * HH * HH));
      const PetscScalar alpha = (g / (2.0 * HH * ice->k)) - 2.0 * HH * HH * beta;
      for (PetscInt k = 0; k < ks; k++) {
        const PetscScalar depth = HH - grid.zlevels[k];
        const PetscScalar Tpmp = ice->triple_point_temp - ice->beta_CC_grad * depth;
        const PetscScalar d2 = depth * depth;

        T[k] = PetscMin(Tpmp,artm(i,j) + alpha * d2 + beta * d2 * d2);

      }
      for (PetscInt k = ks; k < grid.Mz; k++) // above ice
        T[k] = artm(i,j);
      
      if (!do_cold) {
	for (PetscInt k = 0; k < grid.Mz; ++k) {
	  const PetscScalar depth = HH - grid.zlevels[k];
	  const PetscScalar pressure = 
	    EC->getPressureFromDepth(depth);
	  // reuse T to store enthalpy; assume that the ice is cold
	  ierr = EC->getEnthPermissive(T[k], 0.0, pressure, T[k]); CHKERRQ(ierr);
	}
      }

      ierr = result->setInternalColumn(i,j,T); CHKERRQ(ierr);
      
    }
  }
  ierr =     vH.end_access(); CHKERRQ(ierr);
  ierr =   vbed.end_access(); CHKERRQ(ierr);
  ierr =   vGhf.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr =   artm.end_access(); CHKERRQ(ierr);

  delete [] T;

  ierr = result->beginGhostComm(); CHKERRQ(ierr);
  ierr = result->endGhostComm(); CHKERRQ(ierr);

  if (do_cold) {
    ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);
  }

  return 0;
}

