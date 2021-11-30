// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#include "Anomaly.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/io_helpers.hh"

namespace pism {
namespace atmosphere {

Anomaly::Anomaly(IceGrid::ConstPtr g, AtmosphereModel* in)
  : PGivenClimate<PAModifier,AtmosphereModel>(g, in) {
  m_option_prefix  = "-atmosphere_anomaly";

  process_options();

  m_modify_precip = not options::Bool(m_option_prefix + "_temp_only",
                                      "Modify air temperature only,"
                                      " leaving precipitation unchanged");

  // will be de-allocated by the parent's destructor
  {
    m_air_temp_anomaly = new IceModelVec2T;
    m_fields["delta_T"] = m_air_temp_anomaly;
  }

  if (m_modify_precip) {
    m_precipitation_anomaly = new IceModelVec2T;
    m_fields["precipitation_anomaly"] = m_precipitation_anomaly;
  }

  // calls set_n_records(), so ...->create() has to be called after set_vec_parameters()
  set_vec_parameters({});

  {
    m_air_temp_anomaly->create(m_grid, "delta_T");
    m_air_temp_anomaly->set_attrs("climate_forcing",
                                  "anomaly of the near-surface air temperature",
                                  "Kelvin", "");
  }

  if (m_modify_precip) {
    m_precipitation_anomaly->create(m_grid, "precipitation_anomaly");
    m_precipitation_anomaly->set_attrs("climate_forcing",
                                       "anomaly of the ice-equivalent precipitation rate",
                                       "kg m-2 second-1", "");
    m_precipitation_anomaly->metadata().set_string("glaciological_units", "kg m-2 year-1");
  }
}

Anomaly::~Anomaly() {
  // empty
}

void Anomaly::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2,
             "* Initializing the -atmosphere ...,anomaly code...\n"
             "    reading anomalies from %s ...\n",
             m_filename.c_str());

  m_air_temp_anomaly->init(m_filename, m_bc_period, m_bc_reference_time);

  if (m_modify_precip) {
    m_precipitation_anomaly->init(m_filename, m_bc_period, m_bc_reference_time);
  }
}

void Anomaly::update_impl(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  m_air_temp_anomaly->average(m_t, m_dt);

  if (m_modify_precip) {
    m_precipitation_anomaly->average(m_t, m_dt);
  }
}

void Anomaly::mean_precipitation_impl(IceModelVec2S &result) const {
  m_input_model->mean_precipitation(result);

  if (m_modify_precip) {
    result.add(1.0, *m_precipitation_anomaly);
  }
}

void Anomaly::mean_annual_temp_impl(IceModelVec2S &result) const {
  m_input_model->mean_annual_temp(result);

  result.add(1.0, *m_air_temp_anomaly);
}

void Anomaly::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();
  m_air_temp_anomaly->begin_access();

  if (m_modify_precip) {
    m_precipitation_anomaly->begin_access();
  }
}

void Anomaly::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();
  m_air_temp_anomaly->end_access();

  if (m_modify_precip) {
    m_precipitation_anomaly->end_access();
  }
}

void Anomaly::init_timeseries_impl(const std::vector<double> &ts) const {
  PAModifier::init_timeseries_impl(ts);

  m_air_temp_anomaly->init_interpolation(ts);

  if (m_modify_precip) {
    m_precipitation_anomaly->init_interpolation(ts);
  }
}

void Anomaly::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->temp_time_series(i, j, result);

  m_temp_anomaly.reserve(m_ts_times.size());
  m_air_temp_anomaly->interp(i, j, m_temp_anomaly);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_temp_anomaly[k];
  }
}

void Anomaly::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->precip_time_series(i, j, result);

  if (m_modify_precip) {
    m_mass_flux_anomaly.reserve(m_ts_times.size());
    m_precipitation_anomaly->interp(i, j, m_mass_flux_anomaly);

    for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
      result[k] += m_mass_flux_anomaly[k];
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
