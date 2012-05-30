// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir, Andy Aschwanden and Torsten Albrecht
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

#include "POConstantPIK.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

POConstantPIK::POConstantPIK(IceGrid &g, const NCConfigVariable &conf)
  : PISMOceanModel(g, conf) {

  shelfbmassflux.init_2d("shelfbmassflux", g);
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("m s-1");
  shelfbmassflux.set_glaciological_units("m year-1");

  shelfbtemp.init_2d("shelfbtemp", g);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");
}

PetscErrorCode POConstantPIK::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (!config.get_flag("is_dry_simulation")) {
    ierr = verbPrintf(2, grid.com, "* Initializing the constant ocean model...\n"); CHKERRQ(ierr);
  }

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  return 0;
}

PetscErrorCode POConstantPIK::sea_level_elevation(PetscReal &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POConstantPIK::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar T0 = config.get("water_melting_point_temperature"), // K
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME issue #15
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes mass flux in ice-equivalent m s-1.
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
PetscErrorCode POConstantPIK::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  PetscReal L = config.get("water_latent_heat_fusion"),
    rho_ocean = config.get("sea_water_density"),
    rho_ice = config.get("ice_density");

  const PetscScalar c_p_ocean	  = 3974.0,   // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_T	  = 1e-4;     // m/s, thermal exchange velocity
  //FIXME: gamma_T should be a function of the friction velocity, not a const

  PetscScalar ocean_salinity = 35.0; 

  PetscScalar T_water = -1.7, //Default in PISM-PIK
    T_ocean = 273.15 + T_water;

  // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
  // PetscReal meltrate = config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * rho_ice); // m s-1

  PetscReal meltfactor = 5e-3;
  bool meltfactorSet;
  double meltfactor_pik;
  ierr = PISMOptionsReal("-meltfactor_pik",
                         "Uses as a meltfactor as in sub-shelf-melting parameterization of martin_winkelmann11",
                         meltfactor_pik, meltfactorSet); CHKERRQ(ierr);

  if (meltfactorSet) {
    meltfactor = meltfactor_pik; // default is 5e-3 as in martin_winkelmann11
  }

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // compute T_f[i][j] according to beckmann_goosse03, which has the
      // meaning of the freezing temperature of the ocean water directly
      // under the shelf, (of salinity 35psu) [this is related to the
      // Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for
      // details]
      PetscScalar shelfbaseelev = - (rho_ice / rho_ocean) * H[i][j],
        T_f= 273.15 + (0.0939 -0.057 * ocean_salinity + 7.64e-4 * shelfbaseelev);
      // add 273.15 to get it in Kelvin

      // compute ocean_heat_flux according to beckmann_goosse03
      // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
      PetscScalar oceanheatflux = meltfactor * rho_ocean * c_p_ocean * gamma_T * (T_ocean - T_f); // in W/m^2
      // TODO: T_ocean -> field!

      // shelfbmassflux is positive if ice is freezing on; here it is always negative:
      // same sign as OceanHeatFlux... positive if massflux FROM ice TO ocean
      result(i,j) = oceanheatflux / (L * rho_ice); // m s-1

    }
  }

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

void POConstantPIK::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (keyword == "medium" || keyword == "big") {
    result["shelfbtemp"] = shelfbtemp;
    result["shelfbmassflux"] = shelfbmassflux;
  }
}

PetscErrorCode POConstantPIK::define_variables(set<string> vars, const PIO &nc,
                                               PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.define(nc, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode POConstantPIK::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}
