/* Copyright (C) 2016, 2017, 2018 PISM Authors
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
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/coupler/util/init_step.hh"

namespace pism {
namespace ocean {

InitializationHelper::InitializationHelper(IceGrid::ConstPtr g, OceanModel* in)
  : OceanModel(g, in),
    m_sea_level_metadata("effective_sea_level_elevation",
                         m_config->get_string("time.dimension_name"),
                         m_sys) {

  m_melange_back_pressure_fraction = allocate_melange_back_pressure(g);
  m_melange_back_pressure_fraction->set_name("effective_melange_back_pressure_fraction");
  m_melange_back_pressure_fraction->metadata().set_string("pism_intent", "model_state");

  m_shelf_base_temperature = allocate_shelf_base_temperature(g);
  m_shelf_base_temperature->set_name("effective_shelf_base_temperature");
  m_shelf_base_temperature->metadata().set_string("pism_intent", "model_state");

  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);
  m_shelf_base_mass_flux->set_name("effective_shelf_base_mass_flux");
  m_shelf_base_mass_flux->metadata().set_string("pism_intent", "model_state");

  m_sea_level = 0.0;

  m_sea_level_metadata.set_string("pism_intent", "model_state");
  m_sea_level_metadata.set_string("units", "meters");
  m_sea_level_metadata.set_string("long_name", "effective sea level elevation, "
                                  "as seen by the ice dynamics code (for re-starting)");
}

void InitializationHelper::update_impl(double t, double dt) {
  OceanModel::update_impl(t, dt);
}

void InitializationHelper::init_impl() {
  m_input_model->init();

  InputOptions opts = process_input_options(m_grid->com);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "* Reading effective ocean model outputs from '%s' for re-starting...\n",
                   opts.filename.c_str());

    PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
    const unsigned int time_length = file.inq_nrecords();
    const unsigned int last_record = time_length > 0 ? time_length - 1 : 0;

    m_melange_back_pressure_fraction->read(file, last_record);
    m_shelf_base_mass_flux->read(file, last_record);
    m_shelf_base_temperature->read(file, last_record);
    {
      std::vector<double> data;
      file.get_1d_var(m_sea_level_metadata.get_name(),
                      last_record, 1, // start, count
                      data);
      m_sea_level = data[0];
    }

    file.close();
  } else {
    m_log->message(2, "* Performing a 'fake' ocean model time-step for bootstrapping...\n");

    init_step(this, *m_grid->ctx()->time());
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is equivalent to
  // "-i ... -bootstrap -regrid_file ..."
  {
    regrid("ocean model initialization helper", *m_melange_back_pressure_fraction,
           REGRID_WITHOUT_REGRID_VARS);
    regrid("ocean model initialization helper", *m_shelf_base_mass_flux,
           REGRID_WITHOUT_REGRID_VARS);
    regrid("ocean model initialization helper", *m_shelf_base_temperature,
           REGRID_WITHOUT_REGRID_VARS);
  }
  // FIXME: fake "regridding" of sea level
}

void InitializationHelper::define_model_state_impl(const PIO &output) const {
  m_melange_back_pressure_fraction->define(output);
  m_shelf_base_mass_flux->define(output);
  m_shelf_base_temperature->define(output);

  io::define_timeseries(m_sea_level_metadata, output, PISM_DOUBLE);

  m_input_model->define_model_state(output);
}

void InitializationHelper::write_model_state_impl(const PIO &output) const {
  m_melange_back_pressure_fraction->write(output);
  m_shelf_base_mass_flux->write(output);
  m_shelf_base_temperature->write(output);

  const unsigned int
    time_length = output.inq_dimlen(m_sea_level_metadata.get_dimension_name()),
    t_start = time_length > 0 ? time_length - 1 : 0;
  io::write_timeseries(output, m_sea_level_metadata, t_start, m_sea_level,
                       PISM_DOUBLE);

  m_input_model->write_model_state(output);
}

double InitializationHelper::sea_level_elevation_impl() const {
  return m_sea_level;
}

const IceModelVec2S& InitializationHelper::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const IceModelVec2S& InitializationHelper::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

const IceModelVec2S& InitializationHelper::melange_back_pressure_fraction_impl() const {
  return *m_melange_back_pressure_fraction;
}

} // end of namespace ocean
} // end of namespace pism
