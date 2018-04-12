/* Copyright (C) 2018 PISM Authors
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

#include "pism/coupler/util/init_step.hh"

namespace pism {
namespace ocean {
namespace sea_level {

InitializationHelper::InitializationHelper(IceGrid::ConstPtr grid,
                                           std::shared_ptr<SeaLevel> in)
  : SeaLevel(grid, in) {

  m_sea_level.metadata().set_name("effective_sea_level_elevation");
  m_sea_level.metadata().set_string("pism_intent", "model_state");
}

void InitializationHelper::update_impl(const Geometry &geometry, double t, double dt) {
  SeaLevel::update_impl(geometry, t, dt);

  m_sea_level.copy_from(m_input_model->elevation());
}

void InitializationHelper::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "* Reading effective sea level forcing from '%s' for re-starting...\n",
                   opts.filename.c_str());

    PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
    const unsigned int time_length = file.inq_nrecords();
    const unsigned int last_record = time_length > 0 ? time_length - 1 : 0;

    m_sea_level.read(file, last_record);

    file.close();
  } else {
    m_log->message(2, "* Performing a 'fake' sea level forcing time-step for bootstrapping...\n");

    init_step(this, geometry, *m_grid->ctx()->time());
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is
  // equivalent to "-i ... -bootstrap -regrid_file ..."
  {
    regrid("ocean model initialization helper", m_sea_level,
           REGRID_WITHOUT_REGRID_VARS);
  }
}

void InitializationHelper::define_model_state_impl(const PIO &output) const {
  m_sea_level.define(output);

  m_input_model->define_model_state(output);
}

void InitializationHelper::write_model_state_impl(const PIO &output) const {
  m_sea_level.write(output);

  m_input_model->write_model_state(output);
}

const IceModelVec2S& InitializationHelper::sea_level_elevation_impl() const {
  return m_sea_level;
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
