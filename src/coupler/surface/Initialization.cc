/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2023, 2024 PISM Authors
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

#include "pism/coupler/surface/Initialization.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/coupler/util/init_step.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace surface {

InitializationHelper::InitializationHelper(std::shared_ptr<const Grid> grid, std::shared_ptr<SurfaceModel> input)
  : SurfaceModel(grid, input),
    m_mass_flux(m_grid, "effective_climatic_mass_balance"),
    m_temperature(m_grid, "effective_ice_surface_temp")
{

  if (not input) {
    throw RuntimeError(PISM_ERROR_LOCATION, "pism::surface::InitializationHelper got a NULL input model");
  }

  // allocate storage
  {
    m_mass_flux.metadata(0)
        .long_name(
            "surface mass balance (accumulation/ablation) rate, as seen by the ice dynamics code (used for restarting)")
        .units("kg m^-2 s^-1")
        .set_time_independent(false);

    m_temperature.metadata(0)
        .long_name(
            "temperature of the ice at the ice surface but below firn processes, as seen by the ice dynamics code (used for restarting)")
        .units("kelvin")
        .set_time_independent(false);

    m_liquid_water_fraction = allocate_liquid_water_fraction(grid);
    m_liquid_water_fraction->metadata().set_name("effective_ice_surface_liquid_water_fraction");
    m_liquid_water_fraction->metadata(0)
        .long_name(
            "liquid water fraction of the ice at the top surface, as seen by the ice dynamics code (used for restarting)")
        .units("1")
        .set_time_independent(false);

    m_layer_mass = allocate_layer_mass(grid);
    m_layer_mass->metadata().set_name("effective_surface_layer_mass");
    m_layer_mass->metadata(0)
        .long_name(
            "mass held in the surface layer, as seen by the ice dynamics code (used for restarting)")
        .units("kg")
        .set_time_independent(false);

    m_layer_thickness = allocate_layer_thickness(grid);
    m_layer_thickness->metadata().set_name("effective_surface_layer_thickness");
    m_layer_thickness->metadata(0)
        .long_name(
            "thickness of the surface layer, as seen by the ice dynamics code (used for restarting)")
        .units("meters")
        .set_time_independent(false);
  }

  m_accumulation = allocate_accumulation(grid);
  m_accumulation->set_name("effective_" + m_accumulation->get_name());

  m_melt = allocate_melt(grid);
  m_melt->set_name("effective_" + m_melt->get_name());

  m_runoff = allocate_runoff(grid);
  m_runoff->set_name("effective_" + m_runoff->get_name());

  // collect pointers
  m_variables = {&m_mass_flux,
                 &m_temperature,
                 m_liquid_water_fraction.get(),
                 m_layer_mass.get(),
                 m_layer_thickness.get(),
                 m_accumulation.get(),
                 m_melt.get(),
                 m_runoff.get()};
}

void InitializationHelper::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "* Reading effective surface model outputs from '%s' for re-starting...\n",
                   opts.filename.c_str());

    File file(m_grid->com, opts.filename, io::PISM_GUESS, io::PISM_READONLY);
    const unsigned int last_record = file.nrecords() - 1;
    for (auto *v : m_variables) {
      v->read(file, last_record);
    }
  } else {
    m_log->message(2, "* Performing a 'fake' surface model time-step for bootstrapping...\n");

    init_step(this, geometry, time());
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is equivalent to
  // "-i ... -bootstrap -regrid_file ..."
  for (auto *v : m_variables) {
    regrid("surface model initialization helper", *v, REGRID_WITHOUT_REGRID_VARS);
  }
}

void InitializationHelper::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  // store outputs of the input model
  m_mass_flux.copy_from(m_input_model->mass_flux());
  m_temperature.copy_from(m_input_model->temperature());
  m_liquid_water_fraction->copy_from(m_input_model->liquid_water_fraction());
  m_layer_mass->copy_from(m_input_model->layer_mass());
  m_layer_thickness->copy_from(m_input_model->layer_thickness());
  m_accumulation->copy_from(m_input_model->accumulation());
  m_melt->copy_from(m_input_model->melt());
  m_runoff->copy_from(m_input_model->runoff());
}

const array::Scalar &InitializationHelper::layer_thickness_impl() const {
  return *m_layer_thickness;
}

const array::Scalar &InitializationHelper::mass_flux_impl() const {
  return m_mass_flux;
}

const array::Scalar &InitializationHelper::temperature_impl() const {
  return m_temperature;
}

const array::Scalar &InitializationHelper::liquid_water_fraction_impl() const {
  return *m_liquid_water_fraction;
}

const array::Scalar &InitializationHelper::layer_mass_impl() const {
  return *m_layer_mass;
}

const array::Scalar &InitializationHelper::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &InitializationHelper::melt_impl() const {
  return *m_melt;
}

const array::Scalar &InitializationHelper::runoff_impl() const {
  return *m_runoff;
}

void InitializationHelper::define_model_state_impl(const File &output) const {
  for (auto *v : m_variables) {
    v->define(output, io::PISM_DOUBLE);
  }
  m_input_model->define_model_state(output);
}

void InitializationHelper::write_model_state_impl(const File &output) const {
  for (auto *v : m_variables) {
    v->write(output);
  }
  m_input_model->write_model_state(output);
}



} // end of namespace surface
} // end of namespace pism
