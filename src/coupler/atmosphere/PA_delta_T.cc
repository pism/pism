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

#include "PA_delta_T.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace atmosphere {

/// delta_T forcing of near-surface air temperatures

Delta_T::Delta_T(IceGrid::ConstPtr g, AtmosphereModel* in)
  : PScalarForcing<AtmosphereModel,PAModifier>(g, in) {
  m_option_prefix = "-atmosphere_delta_T";
  m_offset_name   = "delta_T";

  m_offset = new Timeseries(*m_grid, m_offset_name, m_config->get_string("time.dimension_name"));
  m_offset->metadata().set_string("units", "Kelvin");
  m_offset->metadata().set_string("long_name", "near-surface air temperature offsets");
  m_offset->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
}

void Delta_T::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2, "* Initializing near-surface air temperature forcing using scalar offsets...\n");

  init_internal();
}

MaxTimestep Delta_T::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere delta_T");
}

void Delta_T::init_timeseries_impl(const std::vector<double> &ts) const {
  PAModifier::init_timeseries_impl(ts);

  m_offset_values.resize(m_ts_times.size());
  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    m_offset_values[k] = (*m_offset)(m_ts_times[k]);
  }
}

void Delta_T::mean_annual_temp_impl(IceModelVec2S &result) const {

  m_input_model->mean_annual_temp(result);
  offset_data(result);
}

void Delta_T::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->temp_time_series(i, j, result);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_offset_values[k];
  }
}

} // end of namespace atmosphere
} // end of namespace pism
