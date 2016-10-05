/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#include "PO_delta_MBP.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace ocean {

Delta_MBP::Delta_MBP(IceGrid::ConstPtr g, OceanModel* in)
  : PScalarForcing<OceanModel,OceanModifier>(g, in) {

  m_option_prefix = "-ocean_delta_MBP";
  m_offset_name   = "delta_MBP";

  m_offset = new Timeseries(*m_grid, m_offset_name, m_config->get_string("time.dimension_name"));

  m_offset->metadata().set_string("units", "1");
  m_offset->metadata().set_string("long_name", "melange back pressure fraction");
  m_offset->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
}

Delta_MBP::~Delta_MBP() {
  // empty
}

void Delta_MBP::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2, "* Initializing melange back pressure fraction forcing...\n");

  init_internal();
}

MaxTimestep Delta_MBP::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Delta_MBP::melange_back_pressure_fraction_impl(IceModelVec2S &result) const {
  m_input_model->melange_back_pressure_fraction(result);

  offset_data(result);
}

} // end of namespace ocean
} // end of namespace pism
