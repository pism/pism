// Copyright (C) 2008-2014 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir, Andy Aschwanden and Torsten Albrecht
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

#include "POConstantPIK.hh"
#include "PISMVars.hh"
#include "PISMConfig.hh"
#include "IceGrid.hh"
#include "iceModelVec.hh"
#include "pism_options.hh"
#include <stdexcept>

namespace pism {

POConstantPIK::POConstantPIK(IceGrid &g, const Config &conf)
  : OceanModel(g, conf),
    shelfbmassflux(g.config.get_unit_system(), "shelfbmassflux", grid),
    shelfbtemp(g.config.get_unit_system(), "shelfbtemp", grid)
{
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("kg m-2 s-1");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");

  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");

  meltfactor = config.get("ocean_pik_melt_factor");
}

POConstantPIK::~POConstantPIK() {
  // empty
}

void POConstantPIK::init(Vars &vars) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, grid.com,
             "* Initializing the constant (PIK) ocean model...\n");

  ice_thickness = vars.get_2d_scalar("land_ice_thickness");

  double meltfactor_pik = meltfactor;
  bool meltfactorSet = false;

  OptionsReal("-meltfactor_pik",
              "Use as a melt factor as in sub-shelf-melting parameterization of [@ref Martinetal2011]",
              meltfactor_pik, meltfactorSet);

  if (meltfactorSet) {
    meltfactor = meltfactor_pik;
  }
}

void POConstantPIK::update(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
}

void POConstantPIK::sea_level_elevation(double &result) {
  result = sea_level;
}

void POConstantPIK::shelf_base_temperature(IceModelVec2S &result) {
  const double
    T0          = config.get("water_melting_point_temperature"), // K
    beta_CC     = config.get("beta_CC"),
    g           = config.get("standard_gravity"),
    ice_density = config.get("ice_density");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(result);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * (*ice_thickness)(i,j); // FIXME task #7297
    // temp is set to melting point at depth
    result(i,j) = T0 - beta_CC * pressure;
  }
}

//! \brief Computes mass flux in [kg m-2 s-1].
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
void POConstantPIK::shelf_base_mass_flux(IceModelVec2S &result) {
  const double
    L                 = config.get("water_latent_heat_fusion"),
    sea_water_density = config.get("sea_water_density"),
    ice_density       = config.get("ice_density"),
    c_p_ocean         = 3974.0, // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_T           = 1e-4,   // m/s, thermal exchange velocity
    ocean_salinity    = 35.0,   // g/kg
    T_ocean           = grid.convert(-1.7, "Celsius", "Kelvin");   //Default in PISM-PIK

  //FIXME: gamma_T should be a function of the friction velocity, not a const

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(result);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // compute T_f[i][j] according to beckmann_goosse03, which has the
    // meaning of the freezing temperature of the ocean water directly
    // under the shelf, (of salinity 35psu) [this is related to the
    // Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for
    // details]
    double
      shelfbaseelev = - (ice_density / sea_water_density) * (*ice_thickness)(i,j),
      T_f           = 273.15 + (0.0939 -0.057 * ocean_salinity + 7.64e-4 * shelfbaseelev);
    // add 273.15 to convert from Celsius to Kelvin

    // compute ocean_heat_flux according to beckmann_goosse03
    // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
    double ocean_heat_flux = meltfactor * sea_water_density * c_p_ocean * gamma_T * (T_ocean - T_f); // in W/m^2
    
    // TODO: T_ocean -> field!

    // shelfbmassflux is positive if ice is freezing on; here it is always negative:
    // same sign as ocean_heat_flux (positive if massflux FROM ice TO ocean)
    result(i,j) = ocean_heat_flux / (L * ice_density); // m s-1

    // convert from [m s-1] to [kg m-2 s-1]:
    result(i,j) *= ice_density;
  }
}

void POConstantPIK::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert("shelfbtemp");
    result.insert("shelfbmassflux");
  }
}

void POConstantPIK::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                               IO_Type nctype) {
  if (set_contains(vars, "shelfbtemp")) {
    shelfbtemp.define(nc, nctype, true);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    shelfbmassflux.define(nc, nctype, true);
  }
}

void POConstantPIK::write_variables(const std::set<std::string> &vars, const PIO &nc) {
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
