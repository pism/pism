/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <gsl/gsl_math.h>

#include "Frac_MBP.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace ocean {

Frac_MBP::Frac_MBP(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in)
  : PScalarForcing<OceanModel,OceanModel>(g, in) {

  m_option_prefix = "-ocean_frac_MBP";
  m_offset_name   = "frac_MBP";

  m_offset.reset(new Timeseries(*m_grid, m_offset_name, m_config->get_string("time.dimension_name")));

  m_offset->variable().set_string("units", "1");
  m_offset->variable().set_string("long_name", "melange back pressure fraction");
  m_offset->dimension().set_string("units", m_grid->ctx()->time()->units_string());

  m_melange_back_pressure_fraction = allocate_melange_back_pressure(g);
}

Frac_MBP::~Frac_MBP() {
  // empty
}

void Frac_MBP::init_impl() {

  m_input_model->init();

  m_log->message(2, "* Initializing melange back pressure fraction forcing...\n");

  init_internal();
}

void Frac_MBP::update_impl(double t, double dt) {
  super::update_impl(t, dt);

  m_melange_back_pressure_fraction->copy_from(m_input_model->melange_back_pressure_fraction());
  m_melange_back_pressure_fraction->scale(m_current_forcing);
}

const IceModelVec2S& Frac_MBP::melange_back_pressure_fraction_impl() const {
  return *m_melange_back_pressure_fraction;
}


} // end of namespace ocean
} // end of namespace pism
