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
namespace frontalmelt {

InitializationHelper::InitializationHelper(IceGrid::ConstPtr g, std::shared_ptr<FrontalMeltModel> in)
  : FrontalMeltModel(g, in) {


  m_frontal_melt_rate = allocate_frontal_melt_rate(g);
  m_frontal_melt_rate->set_name("effective_frontal_melt_rate");
  m_frontal_melt_rate->metadata().set_string("pism_intent", "model_state");
}

void InitializationHelper::update_impl(const Geometry &geometry, double t, double dt) {
  FrontalMeltModel::update_impl(geometry, t, dt);

  m_frontal_melt_rate->copy_from(m_input_model->frontal_melt_rate());
}

void InitializationHelper::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "* Reading effective ocean model outputs from '%s' for re-starting...\n",
                   opts.filename.c_str());

    PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
    const unsigned int time_length = file.inq_nrecords();
    const unsigned int last_record = time_length > 0 ? time_length - 1 : 0;

    m_frontal_melt_rate->read(file, last_record);

    file.close();
  } else {
    m_log->message(2, "* Performing a 'fake' frontal melt model time-step for bootstrapping...\n");

    init_step(this, geometry, *m_grid->ctx()->time());
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is equivalent to
  // "-i ... -bootstrap -regrid_file ..."
  {
    regrid("frontal melt model initialization helper", *m_frontal_melt_rate,
           REGRID_WITHOUT_REGRID_VARS);
  }
}

void InitializationHelper::define_model_state_impl(const PIO &output) const {
  m_frontal_melt_rate->define(output);

  m_input_model->define_model_state(output);
}

void InitializationHelper::write_model_state_impl(const PIO &output) const {
  m_frontal_melt_rate->write(output);

  m_input_model->write_model_state(output);
}

const IceModelVec2S& InitializationHelper::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

} // end of namespace ocean
} // end of namespace pism
