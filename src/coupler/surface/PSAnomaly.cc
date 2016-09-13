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

#include <gsl/gsl_math.h>

#include "PSAnomaly.hh"
#include "base/util/IceGrid.hh"
#include "base/util/io/io_helpers.hh"

namespace pism {
namespace surface {

Anomaly::Anomaly(IceGrid::ConstPtr g, SurfaceModel* in)
  : PGivenClimate<SurfaceModifier,SurfaceModel>(g, in),
    m_climatic_mass_balance(m_sys, "climatic_mass_balance"),
    m_ice_surface_temp(m_sys, "ice_surface_temp") {

  m_option_prefix  = "-surface_anomaly";

  // will be de-allocated by the parent's destructor
  m_climatic_mass_balance_anomaly = new IceModelVec2T;
  m_ice_surface_temp_anomaly      = new IceModelVec2T;

  m_fields["climatic_mass_balance_anomaly"] = m_climatic_mass_balance_anomaly;
  m_fields["ice_surface_temp_anomaly"] = m_ice_surface_temp_anomaly;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  m_ice_surface_temp_anomaly->create(m_grid, "ice_surface_temp_anomaly");
  m_climatic_mass_balance_anomaly->create(m_grid, "climatic_mass_balance_anomaly");

  m_ice_surface_temp_anomaly->set_attrs("climate_forcing",
                                      "anomaly of the temperature of the ice at the ice surface"
                                      " but below firn processes",
                                      "Kelvin", "");
  m_climatic_mass_balance_anomaly->set_attrs("climate_forcing",
                                           "anomaly of the surface mass balance (accumulation/ablation) rate",
                                           "kg m-2 s-1", "");
  m_climatic_mass_balance_anomaly->metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_climatic_mass_balance_anomaly->write_in_glaciological_units = true;

  m_climatic_mass_balance.set_string("pism_intent", "diagnostic");
  m_climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  m_climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.set_string("units", "kg m-2 s-1");
  m_climatic_mass_balance.set_string("glaciological_units", "kg m-2 year-1");

  m_ice_surface_temp.set_string("pism_intent", "diagnostic");
  m_ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  m_ice_surface_temp.set_string("units", "K");
}

Anomaly::~Anomaly() {
  // empty
}

void Anomaly::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (m_input_model != NULL) {
    m_input_model->init();
  }

  m_log->message(2,
             "* Initializing the '-surface ...,anomaly' modifier...\n");

  m_log->message(2,
             "    reading anomalies from %s ...\n", m_filename.c_str());

  m_ice_surface_temp_anomaly->init(m_filename, m_bc_period, m_bc_reference_time);
  m_climatic_mass_balance_anomaly->init(m_filename, m_bc_period, m_bc_reference_time);
}

void Anomaly::update_impl(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  m_climatic_mass_balance_anomaly->average(m_t, m_dt);
  m_ice_surface_temp_anomaly->average(m_t, m_dt);
}

void Anomaly::ice_surface_mass_flux_impl(IceModelVec2S &result) {
  m_input_model->ice_surface_mass_flux(result);
  result.add(1.0, *m_climatic_mass_balance_anomaly);
}

void Anomaly::ice_surface_temperature_impl(IceModelVec2S &result) {
  m_input_model->ice_surface_temperature(result);
  result.add(1.0, *m_ice_surface_temp_anomaly);
}

void Anomaly::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  m_input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big" || keyword == "2dbig") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

void Anomaly::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  std::string order = m_config->get_string("output_variable_order");

  if (set_contains(vars, "ice_surface_temp")) {
    io::define_spatial_variable(m_ice_surface_temp, *m_grid, nc, nctype, order, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    io::define_spatial_variable(m_climatic_mass_balance, *m_grid, nc, nctype, order, true);
  }

  m_input_model->define_variables(vars, nc, nctype);
}

void Anomaly::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = m_ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
    tmp.metadata() = m_climatic_mass_balance;

    ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("climatic_mass_balance");
  }

  m_input_model->write_variables(vars, nc);
}

} // end of namespace surface
} // end of namespace pism
