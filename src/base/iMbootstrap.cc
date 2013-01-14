// Copyright (C) 2004-2013 Jed Brown, Nathan Shemonski, Ed Bueler and
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
#include "PIO.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

//! Read file and use heuristics to initialize PISM from typical 2d data available through remote sensing.
/*! 
This procedure is called by the base class when option <tt>-boot_file</tt> is used.

See chapter 4 of the User's Manual.  We read only 2D information from the bootstrap file.
 */
PetscErrorCode IceModel::bootstrapFromFile(string filename) {
  PetscErrorCode  ierr;

  // Bootstrap 2D fields:
  ierr = bootstrap_2d(filename); CHKERRQ(ierr);

  // Regrid 2D fields:
  ierr = regrid(2); CHKERRQ(ierr);

  // Check the consistency of geometry fields:
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); 

  ierr = verbPrintf(2, grid.com,
		    "getting surface B.C. from couplers...\n"); CHKERRQ(ierr);

  // Update couplers (because heuristics in bootstrap_3d() might need boundary
  // conditions provided by couplers):
  if (surface != NULL) {
    PetscReal max_dt = 0;
    bool restrict = false;
    // FIXME: this will break if a surface model requires contiguous update intervals
    ierr = surface->max_timestep(grid.time->start(), max_dt, restrict); CHKERRQ(ierr);

    if (restrict == false)
      max_dt = convert(1, "year", "seconds");

    ierr = surface->update(grid.time->start(), max_dt); CHKERRQ(ierr);
  } else SETERRQ(grid.com, 1, "surface == NULL");

  if (ocean != NULL) {
    PetscReal max_dt = 0;
    bool restrict = false;
    // FIXME: this will break if an ocean model requires contiguous update intervals
    ierr = ocean->max_timestep(grid.time->start(), max_dt, restrict); CHKERRQ(ierr);

    if (restrict == false)
      max_dt = convert(1, "year", "seconds");

    ierr = ocean->update(grid.time->start(), max_dt); CHKERRQ(ierr);
  } else SETERRQ(grid.com, 1, "ocean == NULL");

  ierr = verbPrintf(2, grid.com,
		    "bootstrapping 3D variables...\n"); CHKERRQ(ierr);

  // Fill 3D fields using heuristics:
  ierr = bootstrap_3d(); CHKERRQ(ierr);

  ierr = regrid(3); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "done reading %s; bootstrapping done\n",filename.c_str()); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::bootstrap_2d(string filename) {
  PetscErrorCode ierr;

  PIO nc(grid, "guess_mode");
  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
		    "bootstrapping by PISM default method from file %s\n", filename.c_str()); CHKERRQ(ierr);

  // report on resulting computational box, rescale grid, actually create a
  // local interpolation context
  ierr = verbPrintf(2, grid.com, 
                    "  rescaling computational box for ice from -boot_file file and\n"
                    "    user options to dimensions:\n"
                    "    [%6.2f km, %6.2f km] x [%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
                    (grid.x0 - grid.Lx)/1000.0,
                    (grid.x0 + grid.Lx)/1000.0,
                    (grid.y0 - grid.Ly)/1000.0,
                    (grid.y0 + grid.Ly)/1000.0,
                    grid.Lz); 
  CHKERRQ(ierr);

  string usurf_name;
  bool hExists=false, maskExists=false, usurf_found_by_std_name;
  ierr = nc.inq_var("usurf", "surface_altitude", hExists, usurf_name, usurf_found_by_std_name); CHKERRQ(ierr);
  ierr = nc.inq_var("mask", maskExists); CHKERRQ(ierr);

  string lon_name, lat_name;
  bool lonExists=false, latExists=false, lon_found_by_std_name, lat_found_by_std_name;
  ierr = nc.inq_var("lon", "longitude", lonExists, lon_name, lon_found_by_std_name); CHKERRQ(ierr);
  ierr = nc.inq_var("lat", "latitude",  latExists, lat_name, lat_found_by_std_name); CHKERRQ(ierr);

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
  if (!lonExists) {
    ierr = vLongitude.set_attr("missing_at_bootstrap","true"); CHKERRQ(ierr);
  }

  ierr =  vLatitude.regrid(filename, false); CHKERRQ(ierr);
  if (!latExists) {
    ierr = vLatitude.set_attr("missing_at_bootstrap","true"); CHKERRQ(ierr);
  }

  ierr =         vH.regrid(filename, 
                           config.get("bootstrapping_H_value_no_var")); CHKERRQ(ierr);
  ierr =       vbed.regrid(filename,  
                           config.get("bootstrapping_bed_value_no_var")); CHKERRQ(ierr);
  ierr =       vbmr.regrid(filename,  
                           config.get("bootstrapping_bmelt_value_no_var")); CHKERRQ(ierr);
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
    if (config.get_flag("part_redist")) {
      ierr = vHresidual.set(0.0); CHKERRQ(ierr);
    }
  }

  if (config.get_flag("kill_icebergs")) {
    // will be updated in updateSurfaceElevationAndMask()
    ierr = vIcebergMask.set(ICEBERGMASK_NOT_SET); CHKERRQ(ierr);
  }

  if (config.get_flag("do_eigen_calving")) {
    ierr = strain_rates.set(0.0); CHKERRQ(ierr);
  }

  if (config.get_flag("ssa_dirichlet_bc")) {
    // Do not use Dirichlet B.C. anywhere if bcflag is not present.
    ierr = vBCMask.regrid(filename, 0.0); CHKERRQ(ierr);
    // In the absence of u_ssa_bc and v_ssa_bc in the file the only B.C. that
    // makes sense is the zero Dirichlet B.C.
    ierr = vBCvel.regrid(filename,  0.0); CHKERRQ(ierr);
  }

  // check if Lz is valid
  PetscReal thk_min, thk_max;
  ierr = vH.range(thk_min, thk_max); CHKERRQ(ierr);

  if (thk_max > grid.Lz) {
    PetscPrintf(grid.com,
                "PISM ERROR: Max. ice thickness (%3.3f m) exceeds the height of the computational domain (%3.3f m).\n"
                "            Exiting...\n", thk_max, grid.Lz);
    PISMEnd();
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
      tau3.set(config.get("initial_age_of_ice_years", "years", "seconds"));
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

//! Create a temperature field within ice and bedrock from given surface temperature and geothermal flux maps.
/*!
In bootstrapping we need to guess about the temperature within the ice and
bedrock if surface temperature and geothermal flux maps are given. This rule is
heuristic but seems to work well anyway. Full bootstrapping will start from the
temperature computed by this procedure and then run for a long time (e.g.
\f$10^5\f$ years), with fixed geometry, to get closer to thermomechanically
coupled equilibrium. See the part of the <i>User's Manual</i> on
EISMINT-Greenland.

Consider a horizontal grid point <tt>i,j</tt>. Suppose the surface temperature
\f$T_s\f$ and the geothermal flux \f$g\f$ are given at that grid point. Within
the corresponding column, denote the temperature by \f$T(z)\f$ for some
elevation \f$z\f$ above the base of the ice. (Note ice corresponds to \f$z>0\f$
while bedrock has \f$z<0\f$.) Apply the rule that \f$T(z)=T_s\f$ is \f$z\f$ is
above the top of the ice (at \f$z=H\f$).

Within the ice, set
\f[T(z) = T_s + \alpha (H-z)^2 + \beta (H-z)^4\f]
where \f$\alpha,\beta\f$ are chosen so that
\f[\frac{\partial T}{\partial z}\Big|_{z=0} = - \frac{g}{k_i}\f]
and 
\f[\frac{\partial T}{\partial z}\Big|_{z=H/4} = - \frac{g}{2 k_i}.\f]

The point of the second condition is our observation that, in observed ice, the
rate of decrease in ice temperature with elevation is significantly decreased
at only one quarter of the ice thickness above the base.

The temperature within the ice is not allowed to exceed the pressure-melting
temperature.

Note that the above heuristic rule for ice determines \f$T(0)\f$. Within the
bedrock our rule is that the rate of change with depth is exactly the
geothermal flux:
\f[T(z) = T(0) - \frac{g}{k_r} z.\f]
Note that \f$z\f$ here is negative, so the temperature increases as one goes
down into the bed.

FIXME issue #15
*/
PetscErrorCode IceModel::putTempAtDepth() {
  PetscErrorCode  ierr;

  PetscScalar *T = new PetscScalar[grid.Mz];
  const bool do_cold = config.get_flag("do_cold_ice_methods");
  const PetscScalar ice_k = config.get("ice_thermal_conductivity"),
    melting_point_temp = config.get("water_melting_point_temperature"),
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity");

  if (surface != NULL) {
    ierr = surface->ice_surface_temperature(artm); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1, "PISM ERROR: surface == NULL");
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
      const PetscScalar beta = (4.0/21.0) * (g / (2.0 * ice_k * HH * HH * HH));
      const PetscScalar alpha = (g / (2.0 * HH * ice_k)) - 2.0 * HH * HH * beta;
      for (PetscInt k = 0; k < ks; k++) {
        const PetscScalar depth = HH - grid.zlevels[k];
        const PetscScalar Tpmp = melting_point_temp - beta_CC_grad * depth;
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

  ierr = result->update_ghosts(); CHKERRQ(ierr);

  if (do_cold) {
    ierr = compute_enthalpy_cold(T3, Enth3); CHKERRQ(ierr);
  }

  return 0;
}

