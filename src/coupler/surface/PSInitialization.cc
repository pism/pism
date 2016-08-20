/* Copyright (C) 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PSInitialization.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

Initialization::Initialization(IceGrid::ConstPtr g, SurfaceModel* in)
  : SurfaceModifier(g, in), m_update_called(false) {

  if (in == NULL) {
    throw RuntimeError("pism::surface::Initialization got a NULL input model");
  }

  m_variables.push_back(&m_ice_surface_mass_flux);
  m_variables.push_back(&m_ice_surface_temperature);
  m_variables.push_back(&m_ice_surface_liquid_water_fraction);
  m_variables.push_back(&m_mass_held_in_surface_layer);
  m_variables.push_back(&m_surface_layer_thickness);
}

void Initialization::attach_atmosphere_model_impl(atmosphere::AtmosphereModel *in) {
  m_input_model->attach_atmosphere_model(in);
}

void Initialization::init_impl() {
  m_input_model->init();
}

void Initialization::update_impl(double t, double dt) {
  // update the input model
  SurfaceModifier::update_impl(t, dt);

  // store outputs of the input model
  m_input_model->ice_surface_mass_flux(m_ice_surface_mass_flux);
  m_input_model->ice_surface_temperature(m_ice_surface_temperature);
  m_input_model->ice_surface_liquid_water_fraction(m_ice_surface_liquid_water_fraction);
  m_input_model->mass_held_in_surface_layer(m_mass_held_in_surface_layer);
  m_input_model->surface_layer_thickness(m_surface_layer_thickness);

  m_update_called = true;
}

void Initialization::ice_surface_mass_flux_impl(IceModelVec2S &result) {
  result.copy_from(m_ice_surface_mass_flux);
}

void Initialization::ice_surface_temperature_impl(IceModelVec2S &result) {
  result.copy_from(m_ice_surface_temperature);
}

void Initialization::ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) {
  result.copy_from(m_ice_surface_liquid_water_fraction);
}

void Initialization::mass_held_in_surface_layer_impl(IceModelVec2S &result) {
  result.copy_from(m_mass_held_in_surface_layer);
}

void Initialization::surface_layer_thickness_impl(IceModelVec2S &result) {
  result.copy_from(m_surface_layer_thickness);
}

void Initialization::add_vars_to_output_impl(const std::string &keyword,
                                             std::set<std::string> &result) {
  // add all the variables we keep track of
  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    result.insert(m_variables[k]->get_name());
  }

  m_input_model->add_vars_to_output(keyword, result);
}

static bool in(const std::set<std::string> &S, const IceModelVec *vec) {
  return set_contains(S, vec->get_name());
}

void Initialization::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                           IO_Type nctype) {
  // make a copy of the set of variables so that we can modify it
  std::set<std::string> list = vars;

  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    const IceModelVec *variable = m_variables[k];
    if (in(list, variable)) {
      variable->define(nc, nctype);
      list.erase(variable->get_name());
    }
  }

  m_input_model->define_variables(list, nc, nctype);
}

void Initialization::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  // make a copy of the set of variables so that we can modify it
  std::set<std::string> list = vars;

  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    const IceModelVec *variable = m_variables[k];
    if (in(list, variable)) {
      variable->write(nc);
      list.erase(variable->get_name());
    }
  }

  m_input_model->write_variables(list, nc);
}

} // end of namespace surface
} // end of namespace pism
