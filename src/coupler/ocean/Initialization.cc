/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#include "Initialization.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_options.hh"
#include "pism/coupler/util/init_step.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace ocean {

InitializationHelper::InitializationHelper(std::shared_ptr<const Grid> g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_water_column_pressure = allocate_water_column_pressure(g);
  m_water_column_pressure->set_name("effective_water_column_pressure");
  m_water_column_pressure->metadata()["pism_intent"] = "model_state";

  m_shelf_base_temperature = allocate_shelf_base_temperature(g);
  m_shelf_base_temperature->set_name("effective_shelf_base_temperature");
  m_shelf_base_temperature->metadata()["pism_intent"] = "model_state";

  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);
  m_shelf_base_mass_flux->set_name("effective_shelf_base_mass_flux");
  // use internal units when saving
  auto units = m_shelf_base_mass_flux->metadata().get_string("units");
  m_shelf_base_mass_flux->metadata()["glaciological_units"] = units;
  m_shelf_base_mass_flux->metadata()["pism_intent"] = "model_state";
}

void InitializationHelper::update_impl(const Geometry &geometry, double t, double dt) {
  OceanModel::update_impl(geometry, t, dt);

  m_water_column_pressure->copy_from(m_input_model->average_water_column_pressure());
  m_shelf_base_temperature->copy_from(m_input_model->shelf_base_temperature());
  m_shelf_base_mass_flux->copy_from(m_input_model->shelf_base_mass_flux());
}

void InitializationHelper::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "* Reading effective ocean model outputs from '%s' for re-starting...\n",
                   opts.filename.c_str());

    File file(m_grid->com, opts.filename, io::PISM_GUESS, io::PISM_READONLY);
    const unsigned int time_length = file.nrecords();
    const unsigned int last_record = time_length > 0 ? time_length - 1 : 0;

    m_water_column_pressure->read(file, last_record);
    m_shelf_base_mass_flux->read(file, last_record);
    m_shelf_base_temperature->read(file, last_record);

    file.close();
  } else {
    m_log->message(2, "* Performing a 'fake' ocean model time-step for bootstrapping...\n");

    init_step(this, geometry, time());
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is equivalent to
  // "-i ... -bootstrap -regrid_file ..."
  {
    regrid("ocean model initialization helper", *m_water_column_pressure,
           REGRID_WITHOUT_REGRID_VARS);
    regrid("ocean model initialization helper", *m_shelf_base_mass_flux,
           REGRID_WITHOUT_REGRID_VARS);
    regrid("ocean model initialization helper", *m_shelf_base_temperature,
           REGRID_WITHOUT_REGRID_VARS);
  }
}

void InitializationHelper::define_model_state_impl(const File &output) const {
  m_water_column_pressure->define(output, io::PISM_DOUBLE);
  m_shelf_base_mass_flux->define(output, io::PISM_DOUBLE);
  m_shelf_base_temperature->define(output, io::PISM_DOUBLE);

  m_input_model->define_model_state(output);
}

void InitializationHelper::write_model_state_impl(const File &output) const {
  m_water_column_pressure->write(output);
  m_shelf_base_mass_flux->write(output);
  m_shelf_base_temperature->write(output);

  m_input_model->write_model_state(output);
}

const array::Scalar& InitializationHelper::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const array::Scalar& InitializationHelper::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

const array::Scalar& InitializationHelper::average_water_column_pressure_impl() const {
  return *m_water_column_pressure;
}

} // end of namespace ocean
} // end of namespace pism
