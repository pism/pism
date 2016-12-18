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
  : Component_TS(g) {
  // empty
}

AtmosphereModel::~AtmosphereModel() {
  // empty
}

void AtmosphereModel::init() {
  this->init_impl();
}

void AtmosphereModel::mean_precipitation(IceModelVec2S &result) const {
  this->mean_precipitation_impl(result);
}

void AtmosphereModel::mean_annual_temp(IceModelVec2S &result) const {
  this->mean_annual_temp_impl(result);
}

void AtmosphereModel::begin_pointwise_access() const {
  this->begin_pointwise_access_impl();
}

void AtmosphereModel::end_pointwise_access() const {
  this->end_pointwise_access_impl();
}

void AtmosphereModel::init_timeseries(const std::vector<double> &ts) const {
  this->init_timeseries_impl(ts);
}

void AtmosphereModel::precip_time_series(int i, int j, std::vector<double> &result) const {
  this->precip_time_series_impl(i, j, result);
}

void AtmosphereModel::temp_time_series(int i, int j, std::vector<double> &result) const {
  this->temp_time_series_impl(i, j, result);
}

std::map<std::string, Diagnostic::Ptr> AtmosphereModel::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"air_temp_snapshot",       Diagnostic::Ptr(new PA_air_temp_snapshot(this))},
    {"effective_air_temp",      Diagnostic::Ptr(new PA_air_temp(this))},
    {"effective_precipitation", Diagnostic::Ptr(new PA_precipitation(this))},
  };
  return result;
}

PA_air_temp_snapshot::PA_air_temp_snapshot(const AtmosphereModel *m)
  : Diag<AtmosphereModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "air_temp_snapshot")};

  set_attrs("instantaneous value of the near-surface air temperature",
            "",                 // no standard name
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PA_air_temp_snapshot::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "air_temp_snapshot", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  std::vector<double> current_time(1, m_grid->ctx()->time()->current());
  std::vector<double> temperature(1, 0.0);

  model->init_timeseries(current_time);

  model->begin_pointwise_access();

  IceModelVec::AccessList list(*result);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      model->temp_time_series(i, j, temperature);

      (*result)(i, j) = temperature[0];
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  model->end_pointwise_access();

  return result;
}

PA_air_temp::PA_air_temp(const AtmosphereModel *m)
  : Diag<AtmosphereModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "effective_air_temp")};

  set_attrs("effective mean-annual near-surface air temperature", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PA_air_temp::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "effective_air_temp", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->mean_annual_temp(*result);

  return result;
}

PA_precipitation::PA_precipitation(const AtmosphereModel *m)
  : Diag<AtmosphereModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "effective_precipitation")};

  set_attrs("effective precipitation rate",
            "",                 // no standard name, as far as I know
            "kg m-2 second-1", "kg m-2 year-1", 0);
}

IceModelVec::Ptr PA_precipitation::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "effective_precipitation", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->mean_precipitation(*result);

  return result;
}

} // end of namespace atmosphere
} // end of namespace pism
