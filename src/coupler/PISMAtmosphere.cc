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

#include "PISMAtmosphere.hh"
#include "base/util/error_handling.hh"

namespace pism {
namespace atmosphere {

AtmosphereModel::AtmosphereModel(IceGrid::ConstPtr g)
  : Component_TS(g),
    m_air_temp(m_sys, "effective_air_temp"),
    m_precipitation(m_sys, "effective_precipitation") {

  m_air_temp.set_string("pism_intent", "diagnostic");
  m_air_temp.set_string("long_name", "near-surface air temperature");
  m_air_temp.set_string("units", "K");

  m_precipitation.set_string("pism_intent", "diagnostic");
  m_precipitation.set_string("long_name", "precipitation");
  m_precipitation.set_string("units", "kg m-2 second-1");
  m_precipitation.set_string("glaciological_units", "kg m-2 year-1");
}

AtmosphereModel::~AtmosphereModel() {
  // empty
}

void AtmosphereModel::init() {
  this->init_impl();
}

void AtmosphereModel::mean_precipitation(IceModelVec2S &result) {
  this->mean_precipitation_impl(result);
}

void AtmosphereModel::mean_annual_temp(IceModelVec2S &result) {
  this->mean_annual_temp_impl(result);
}

void AtmosphereModel::begin_pointwise_access() {
  this->begin_pointwise_access_impl();
}

void AtmosphereModel::end_pointwise_access() {
  this->end_pointwise_access_impl();
}

void AtmosphereModel::init_timeseries(const std::vector<double> &ts) {
  this->init_timeseries_impl(ts);
}

void AtmosphereModel::precip_time_series(int i, int j, std::vector<double> &result) {
  this->precip_time_series_impl(i, j, result);
}

void AtmosphereModel::temp_time_series(int i, int j, std::vector<double> &result) {
  this->temp_time_series_impl(i, j, result);
}

void AtmosphereModel::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                           std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
  // Don't override the diagnostic if it is already in dict.
  if (not dict["air_temp_snapshot"]) {
    dict["air_temp_snapshot"] = Diagnostic::Ptr(new PA_air_temp_snapshot(this));
  }
  (void) ts_dict;
}

PA_air_temp_snapshot::PA_air_temp_snapshot(AtmosphereModel *m)
  : Diag<AtmosphereModel>(m) {

  /* set metadata: */
  m_vars.push_back(SpatialVariableMetadata(m_sys, "air_temp_snapshot"));

  set_attrs("instantaneous value of the near-surface air temperature",
            "",                 // no standard name
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PA_air_temp_snapshot::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "air_temp_snapshot", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  std::vector<double> current_time(1, m_grid->ctx()->time()->current());
  std::vector<double> temp(1, 0.0);

  model->init_timeseries(current_time);

  model->begin_pointwise_access();

  IceModelVec::AccessList list(*result);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      model->temp_time_series(i, j, temp);

      (*result)(i, j) = temp[0];
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  model->end_pointwise_access();

  return result;
}

} // end of namespace atmosphere
} // end of namespace pism
