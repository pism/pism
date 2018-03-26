// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "Delta_T.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/coupler/util/ScalarForcing.hh"

namespace pism {
namespace atmosphere {

/// delta_T forcing of near-surface air temperatures

Delta_T::Delta_T(IceGrid::ConstPtr grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in) {

  m_forcing.reset(new ScalarForcing(grid->ctx(),
                                    "-atmosphere_delta_T",
                                    "delta_T",
                                    "Kelvin",
                                    "Kelvin",
                                    "near-surface air temperature offsets"));

  m_temperature = allocate_temperature(grid);
}

void Delta_T::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2, "* Initializing near-surface air temperature forcing using scalar offsets...\n");

  m_forcing->init();
}

void Delta_T::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_offset_values.resize(ts.size());
  for (unsigned int k = 0; k < ts.size(); ++k) {
    m_offset_values[k] = m_forcing->value(ts[k]);
  }
}

void Delta_T::update_impl(const Geometry& geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);
  m_forcing->update(t, dt);

  m_temperature->copy_from(m_input_model->mean_annual_temp());
  m_temperature->shift(m_forcing->value());
}

const IceModelVec2S& Delta_T::mean_annual_temp_impl() const {
  return *m_temperature;
}

void Delta_T::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->temp_time_series(i, j, result);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_offset_values[k];
  }
}

} // end of namespace atmosphere
} // end of namespace pism
