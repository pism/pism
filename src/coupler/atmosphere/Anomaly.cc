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

#include "Anomaly.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/coupler/util/options.hh"

namespace pism {
namespace atmosphere {

Anomaly::Anomaly(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(g, in) {

  ForcingOptions opt(*m_grid->ctx(), "atmosphere.anomaly");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_air_temp_anomaly = IceModelVec2T::ForcingField(m_grid,
                                                     file,
                                                     "air_temp_anomaly",
                                                     "", // no standard name
                                                     buffer_size,
                                                     evaluations_per_year,
                                                     periodic);

    m_precipitation_anomaly = IceModelVec2T::ForcingField(m_grid,
                                                          file,
                                                          "precipitation_anomaly",
                                                          "", // no standard name
                                                          buffer_size,
                                                          evaluations_per_year,
                                                          periodic);
  }

  m_air_temp_anomaly->set_attrs("climate_forcing",
                                "anomaly of the near-surface air temperature",
                                "Kelvin", "");

  m_precipitation_anomaly->set_attrs("climate_forcing",
                                     "anomaly of the ice-equivalent precipitation rate",
                                     "kg m-2 second-1", "");
  m_precipitation_anomaly->metadata().set_string("glaciological_units", "kg m-2 year-1");

  m_precipitation = allocate_precipitation(g);
  m_temperature   = allocate_temperature(g);
}

Anomaly::~Anomaly()
{
  // empty
}

void Anomaly::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  ForcingOptions opt(*m_grid->ctx(), "atmosphere.anomaly");

  m_log->message(2,
                 "* Initializing the -atmosphere ...,anomaly code...\n");

  m_log->message(2,
                 "    reading anomalies from %s ...\n",
                 opt.filename.c_str());

  m_air_temp_anomaly->init(opt.filename, opt.period, opt.reference_time);
  m_precipitation_anomaly->init(opt.filename, opt.period, opt.reference_time);
}

void Anomaly::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_precipitation_anomaly->update(t, dt);
  m_air_temp_anomaly->update(t, dt);

  m_precipitation_anomaly->average(t, dt);
  m_air_temp_anomaly->average(t, dt);

  // precipitation
  {
    m_precipitation->copy_from(m_input_model->mean_precipitation());
    m_precipitation->add(1.0, *m_precipitation_anomaly);
  }

  // temperature
  {
    m_temperature->copy_from(m_input_model->mean_annual_temp());
    m_temperature->add(1.0, *m_air_temp_anomaly);
  }
}

const IceModelVec2S& Anomaly::mean_precipitation_impl() const {
  return *m_precipitation;
}

const IceModelVec2S& Anomaly::mean_annual_temp_impl() const {
  return *m_temperature;
}

void Anomaly::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();
  m_air_temp_anomaly->begin_access();
  m_precipitation_anomaly->begin_access();
}

void Anomaly::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();
  m_precipitation_anomaly->end_access();
  m_air_temp_anomaly->end_access();
}

void Anomaly::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_air_temp_anomaly->init_interpolation(ts);

  m_precipitation_anomaly->init_interpolation(ts);
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

  m_mass_flux_anomaly.reserve(m_ts_times.size());
  m_precipitation_anomaly->interp(i, j, m_mass_flux_anomaly);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_mass_flux_anomaly[k];
  }
}

} // end of namespace atmosphere
} // end of namespace pism
