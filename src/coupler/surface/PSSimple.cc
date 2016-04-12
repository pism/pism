// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <cassert>
#include <gsl/gsl_math.h>

#include "PSSimple.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_const.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

///// Simple PISM surface model.
Simple::Simple(IceGrid::ConstPtr g)
  : SurfaceModel(g),
    climatic_mass_balance(m_sys, "climatic_mass_balance"),
    ice_surface_temp(m_sys, "ice_surface_temp") {

  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_string("units", "kg m-2 s-1");
  climatic_mass_balance.set_string("glaciological_units", "kg m-2 year-1");

  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ice_surface_temp.set_string("units", "K");
}


void Simple::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  assert(m_atmosphere != NULL);
  m_atmosphere->init();

  m_log->message(2,
             "* Initializing the simplest PISM surface (snow) processes model Simple.\n"
             "  It passes atmospheric state directly to upper ice fluid surface:\n"
             "    surface mass balance          := precipitation,\n"
             "    ice upper surface temperature := 2m air temperature.\n");
}

MaxTimestep Simple::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Simple::update_impl(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
  if (m_atmosphere) {
    m_atmosphere->update(my_t, my_dt);
  }
}


void Simple::ice_surface_mass_flux_impl(IceModelVec2S &result) {
  m_atmosphere->mean_precipitation(result);
  result.scale(m_config->get_double("ice_density"));
}

void Simple::ice_surface_temperature_impl(IceModelVec2S &result) {
  m_atmosphere->mean_annual_temp(result);
}

void Simple::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  SurfaceModel::add_vars_to_output_impl(keyword, result);

  if (keyword == "medium" || keyword == "big" || keyword == "2dbig") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

void Simple::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  std::string order = m_config->get_string("output_variable_order");

  if (set_contains(vars, "ice_surface_temp")) {
    io::define_spatial_variable(ice_surface_temp, *m_grid, nc, nctype, order, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    io::define_spatial_variable(climatic_mass_balance, *m_grid, nc, nctype, order, true);
  }

  SurfaceModel::define_variables_impl(vars, nc, nctype);
}

void Simple::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
    tmp.metadata() = climatic_mass_balance;

    ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("climatic_mass_balance");
  }

  SurfaceModel::write_variables_impl(vars, nc);
}

} // end of namespace surface
} // end of namespace pism
