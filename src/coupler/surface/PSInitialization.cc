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
#include "base/util/pism_options.hh"
#include "base/util/io/PIO.hh"

namespace pism {
namespace surface {

InitializationHelper::InitializationHelper(IceGrid::ConstPtr g, SurfaceModel* in)
  : SurfaceModifier(g, in) {

  if (in == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "pism::surface::InitializationHelper got a NULL input model");
  }

  // allocate storage
  {
    m_ice_surface_mass_flux.create(m_grid, "effective_climatic_mass_balance", WITHOUT_GHOSTS);
    m_ice_surface_mass_flux.set_attrs("model_state",
                                      "surface mass balance (accumulation/ablation) rate, as seen by the ice dynamics code (used for restarting)",
                                      "kg m-2 s-1", "");
    m_ice_surface_mass_flux.set_time_independent(false);
    m_ice_surface_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");

    m_ice_surface_temperature.create(m_grid, "effective_ice_surface_temp", WITHOUT_GHOSTS);
    m_ice_surface_temperature.set_attrs("model_state",
                                        "temperature of the ice at the ice surface but below firn processes, as seen by the ice dynamics code (used for restarting)",
                                        "Kelvin", "");
    m_ice_surface_temperature.set_time_independent(false);

    m_ice_surface_liquid_water_fraction.create(m_grid, "effective_ice_surface_liquid_water_fraction", WITHOUT_GHOSTS);
    m_ice_surface_liquid_water_fraction.set_attrs("model_state",
                                                  "liquid water fraction of the ice at the top surface, as seen by the ice dynamics code (used for restarting)",
                                                  "1", "");
    m_ice_surface_liquid_water_fraction.set_time_independent(false);

    m_mass_held_in_surface_layer.create(m_grid, "effective_mass_held_in_surface_layer", WITHOUT_GHOSTS);
    m_mass_held_in_surface_layer.set_attrs("model_state",
                                           "mass held in the surface layer, as seen by the ice dynamics code (used for restarting)",
                                           "kg",
                                           "");
    m_mass_held_in_surface_layer.set_time_independent(false);

    m_surface_layer_thickness.create(m_grid, "effective_surface_layer_thickness", WITHOUT_GHOSTS);
    m_surface_layer_thickness.set_attrs("model_state",
                                        "thickness of the surface layer, as seen by the ice dynamics code (used for restarting)",
                                        "meters", "");
    m_surface_layer_thickness.set_time_independent(false);
  }

  // collect pointers
  m_variables = {&m_ice_surface_mass_flux,
                 &m_ice_surface_temperature,
                 &m_ice_surface_liquid_water_fraction,
                 &m_mass_held_in_surface_layer,
                 &m_surface_layer_thickness};
}

void InitializationHelper::attach_atmosphere_model_impl(atmosphere::AtmosphereModel *in) {
  m_input_model->attach_atmosphere_model(in);
}

void InitializationHelper::init_impl() {
  m_input_model->init();

  InputOptions opts = process_input_options(m_grid->com);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "* Reading effective surface model outputs from '%s' for re-starting...\n",
                   opts.filename.c_str());

    PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
    const unsigned int last_record = file.inq_nrecords() - 1;
    for (auto v : m_variables) {
      v->read(file, last_record);
    }
  } else {
    m_log->message(2, "* Performing a 'fake' surface model time-step for bootstrapping...\n");

    init_step(*this, *m_grid->ctx()->time());
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is equivalent to
  // "-i ... -bootstrap -regrid_file ..."
  for (auto v : m_variables) {
    regrid("surface model initialization helper", *v, REGRID_WITHOUT_REGRID_VARS);
  }
}

void InitializationHelper::update_impl(double t, double dt) {
  // update the input model
  SurfaceModifier::update_impl(t, dt);

  // store outputs of the input model
  m_input_model->ice_surface_mass_flux(m_ice_surface_mass_flux);
  m_input_model->ice_surface_temperature(m_ice_surface_temperature);
  m_input_model->ice_surface_liquid_water_fraction(m_ice_surface_liquid_water_fraction);
  m_input_model->mass_held_in_surface_layer(m_mass_held_in_surface_layer);
  m_input_model->surface_layer_thickness(m_surface_layer_thickness);
}

void InitializationHelper::ice_surface_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_ice_surface_mass_flux);
}

void InitializationHelper::ice_surface_temperature_impl(IceModelVec2S &result) const {
  result.copy_from(m_ice_surface_temperature);
}

void InitializationHelper::ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) const {
  result.copy_from(m_ice_surface_liquid_water_fraction);
}

void InitializationHelper::mass_held_in_surface_layer_impl(IceModelVec2S &result) const {
  result.copy_from(m_mass_held_in_surface_layer);
}

void InitializationHelper::surface_layer_thickness_impl(IceModelVec2S &result) const {
  result.copy_from(m_surface_layer_thickness);
}

void InitializationHelper::define_model_state_impl(const PIO &output) const {
  for (auto v : m_variables) {
    v->define(output);
  }
  m_input_model->define_model_state(output);
}

void InitializationHelper::write_model_state_impl(const PIO &output) const {
  for (auto v : m_variables) {
    v->write(output);
  }
  m_input_model->write_model_state(output);
}

} // end of namespace surface
} // end of namespace pism
