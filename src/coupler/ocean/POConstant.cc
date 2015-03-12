// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#include <gsl/gsl_math.h>

#include "POConstant.hh"
#include "PISMVars.hh"
#include "PISMConfig.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include "error_handling.hh"

namespace pism {
namespace ocean {
Constant::Constant(const IceGrid &g)
  : OceanModel(g),
    m_shelfbmassflux(g.config.get_unit_system(), "shelfbmassflux", m_grid),
    m_shelfbtemp(g.config.get_unit_system(), "shelfbtemp", m_grid) {

  m_mymeltrate = 0.0;
  m_meltrate_set = false;

  m_shelfbmassflux.set_string("pism_intent", "climate_state");
  m_shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  m_shelfbmassflux.set_string("units", "kg m-2 s-1");
  m_shelfbmassflux.set_string("glaciological_units", "kg m-2 year-1");

  m_shelfbtemp.set_string("pism_intent", "climate_state");
  m_shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  m_shelfbtemp.set_string("units", "Kelvin");
}

Constant::~Constant() {
  // empty
}

void Constant::update_impl(double my_t, double my_dt) {
  // do nothing
  m_t = my_t;
  m_dt = my_dt;
}

void Constant::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (!m_config.get_flag("is_dry_simulation")) {
    verbPrintf(2, m_grid.com, "* Initializing the constant ocean model...\n");
  }

  options::Real meltrate("-shelf_base_melt_rate",
                         "Specifies a sub shelf ice-equivalent melt rate in meters/year",
                         m_mymeltrate);

  if (meltrate.is_set()) {
    m_mymeltrate = meltrate;
    verbPrintf(2, m_grid.com,
               "    - option '-shelf_base_melt_rate' seen, "
               "setting basal sub shelf basal melt rate to %.2f m/year ... \n",
               m_mymeltrate);
  }
}

MaxTimestep Constant::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Constant::sea_level_elevation_impl(double &result) {
  result = m_sea_level;
}

void Constant::shelf_base_temperature_impl(IceModelVec2S &result) {
  const double T0 = m_config.get("water_melting_point_temperature"), // K
    beta_CC       = m_config.get("beta_CC"),
    g             = m_config.get("standard_gravity"),
    ice_density   = m_config.get("ice_density");

  const IceModelVec2S *ice_thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(result);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * (*ice_thickness)(i,j); // FIXME issue #15
    // temp is set to melting point at depth
    result(i,j) = T0 - beta_CC * pressure;
  }
}

//! @brief Computes mass flux in [kg m-2 s-1], from assumption that
//! basal heat flux rate converts to mass flux.
void Constant::shelf_base_mass_flux_impl(IceModelVec2S &result) {
  double
    L           = m_config.get("water_latent_heat_fusion"),
    ice_density = m_config.get("ice_density"),
    meltrate    = 0.0;

  if (m_meltrate_set) {

    meltrate = m_grid.convert(m_mymeltrate, "m year-1", "m s-1");

  } else {

    // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
    meltrate = m_config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * ice_density); // m s-1

  }

  // convert to [kg m-2 s-1] = [m s-1] * [kg m-3]
  meltrate = meltrate * ice_density;

  result.set(meltrate);
}

void Constant::add_vars_to_output_impl(const std::string&, std::set<std::string> &result) {
  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

void Constant::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                  IO_Type nctype) {

  if (set_contains(vars, "shelfbtemp")) {
    m_shelfbtemp.define(nc, nctype, true);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    m_shelfbmassflux.define(nc, nctype, true);
  }
}

void Constant::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = m_shelfbtemp;
    shelf_base_temperature(tmp);
    tmp.write(nc);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = m_shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    shelf_base_mass_flux(tmp);
    tmp.write(nc);
  }
}
} // end of namespape ocean
} // end of namespace pism
