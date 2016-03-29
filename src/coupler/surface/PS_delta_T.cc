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

#include "PS_delta_T.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

/// -surface ...,delta_T (scalar forcing of ice surface temperatures)

Delta_T::Delta_T(IceGrid::ConstPtr g, SurfaceModel* in)
  : PScalarForcing<SurfaceModel,SurfaceModifier>(g, in),
    climatic_mass_balance(m_sys, "climatic_mass_balance"),
    ice_surface_temp(m_sys, "ice_surface_temp") {

  option_prefix = "-surface_delta_T";
  offset_name   = "delta_T";

  offset = new Timeseries(*m_grid, offset_name, m_config->get_string("time_dimension_name"));

  offset->metadata().set_string("units", "Kelvin");
  offset->metadata().set_string("long_name", "ice-surface temperature offsets");
  offset->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());

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

Delta_T::~Delta_T() {
  // empty
}

void Delta_T::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  input_model->init();

  m_log->message(2,
             "* Initializing ice-surface temperature forcing using scalar offsets...\n");

  init_internal();
}

MaxTimestep Delta_T::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Delta_T::ice_surface_temperature_impl(IceModelVec2S &result) {
  input_model->ice_surface_temperature(result);
  offset_data(result);
}

void Delta_T::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big" || keyword == "2dbig") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

void Delta_T::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  std::string order = m_config->get_string("output_variable_order");

  if (set_contains(vars, "ice_surface_temp")) {
    io::define_spatial_variable(ice_surface_temp, *m_grid, nc, nctype, order, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    io::define_spatial_variable(climatic_mass_balance, *m_grid, nc, nctype, order, true);
  }

  input_model->define_variables(vars, nc, nctype);
}

void Delta_T::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
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

  input_model->write_variables(vars, nc);
}

} // end of namespace surface
} // end of namespace pism
