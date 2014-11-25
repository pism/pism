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
#include "iceModelVec.hh"
#include <stdexcept>
#include "error_handling.hh"

namespace pism {

POConstant::POConstant(IceGrid &g, const Config &conf)
  : OceanModel(g, conf),
    shelfbmassflux(g.get_unit_system(), "shelfbmassflux", grid),
    shelfbtemp(g.get_unit_system(), "shelfbtemp", grid) {

  mymeltrate = 0.0;
  meltrate_set = false;

  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("kg m-2 s-1");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");

  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");
}

void POConstant::init(Vars &vars) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (!config.get_flag("is_dry_simulation")) {
    verbPrintf(2, grid.com, "* Initializing the constant ocean model...\n");
  }

  OptionsReal("-shelf_base_melt_rate",
              "Specifies a sub shelf ice-equivalent melt rate in meters/year",
              mymeltrate, meltrate_set);

  if (meltrate_set) {
    verbPrintf(2, grid.com,
               "    - option '-shelf_base_melt_rate' seen, "
               "setting basal sub shelf basal melt rate to %.2f m/year ... \n",
               mymeltrate);
  }

  ice_thickness = vars.get_2d_scalar("land_ice_thickness");
}

void POConstant::sea_level_elevation(double &result) {
  result = sea_level;
}

void POConstant::shelf_base_temperature(IceModelVec2S &result) {
  const double T0 = config.get("water_melting_point_temperature"), // K
    beta_CC       = config.get("beta_CC"),
    g             = config.get("standard_gravity"),
    ice_density   = config.get("ice_density");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(result);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * (*ice_thickness)(i,j); // FIXME issue #15
    // temp is set to melting point at depth
    result(i,j) = T0 - beta_CC * pressure;
  }
}

//! @brief Computes mass flux in [kg m-2 s-1], from assumption that
//! basal heat flux rate converts to mass flux.
void POConstant::shelf_base_mass_flux(IceModelVec2S &result) {
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

  result.set(meltrate);
}

void POConstant::add_vars_to_output(const std::string&, std::set<std::string> &result) {
  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

void POConstant::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                  IO_Type nctype) {

  if (set_contains(vars, "shelfbtemp")) {
    shelfbtemp.define(nc, nctype, true);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    shelfbmassflux.define(nc, nctype, true);
  }
}

void POConstant::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      tmp.create(grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = shelfbtemp;
    shelf_base_temperature(tmp);
    tmp.write(nc);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      tmp.create(grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    shelf_base_mass_flux(tmp);
    tmp.write(nc);
  }
}

} // end of namespace pism
