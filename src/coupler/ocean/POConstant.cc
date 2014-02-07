// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "POConstant.hh"
#include "PISMVars.hh"
#include "PISMConfig.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

POConstant::POConstant(IceGrid &g, const PISMConfig &conf)
  : PISMOceanModel(g, conf),
    shelfbmassflux(g.get_unit_system()),
    shelfbtemp(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_POConstant(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PetscErrorCode POConstant::allocate_POConstant() {
  mymeltrate = 0.0;
  meltrate_set = false;

  shelfbmassflux.init_2d("shelfbmassflux", grid);
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("kg m-2 s-1");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");

  shelfbtemp.init_2d("shelfbtemp", grid);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");

  return 0;
}

PetscErrorCode POConstant::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (!config.get_flag("is_dry_simulation")) {
    ierr = verbPrintf(2, grid.com, "* Initializing the constant ocean model...\n"); CHKERRQ(ierr);
  }

  ierr = PetscOptionsBegin(grid.com, "", "Ocean model", ""); CHKERRQ(ierr);

  ierr = PISMOptionsReal("-shelf_base_melt_rate",
                          "Specifies a sub shelf ice-equivalent melt rate in meters/year",
			  mymeltrate, meltrate_set); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (meltrate_set) {
    ierr = verbPrintf(2, grid.com,
                      "    - option '-shelf_base_melt_rate' seen, "
                      "setting basal sub shelf basal melt rate to %.2f m/year ... \n",
                      mymeltrate); CHKERRQ(ierr);
  }

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  return 0;
}

PetscErrorCode POConstant::sea_level_elevation(double &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POConstant::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const double T0 = config.get("water_melting_point_temperature"), // K
    beta_CC       = config.get("beta_CC"),
    g             = config.get("standard_gravity"),
    ice_density   = config.get("ice_density");

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const double pressure = ice_density * g * (*ice_thickness)(i,j); // FIXME issue #15

      // temp is set to melting point at depth
      result(i,j) = T0 - beta_CC * pressure;
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

//! @brief Computes mass flux in [kg m-2 s-1], from assumption that
//! basal heat flux rate converts to mass flux.
PetscErrorCode POConstant::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  double
    L           = config.get("water_latent_heat_fusion"),
    ice_density = config.get("ice_density"),
    meltrate    = 0.0;

  if (meltrate_set) {

    meltrate = grid.convert(mymeltrate, "m year-1", "m s-1");

  } else {

    // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
    meltrate = config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * ice_density); // m s-1

  }

  // convert to [kg m-2 s-1] = [m s-1] * [kg m-3]
  meltrate = meltrate * ice_density;

  ierr = result.set(meltrate); CHKERRQ(ierr);

  return 0;
}

void POConstant::add_vars_to_output(std::string, std::set<std::string> &result) {
  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

PetscErrorCode POConstant::define_variables(std::set<std::string> vars, const PIO &nc,
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

PetscErrorCode POConstant::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    }

    tmp.metadata() = shelfbtemp;
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    }

    tmp.metadata() = shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
  }

  return 0;
}
