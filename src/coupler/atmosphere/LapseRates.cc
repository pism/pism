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

#include "LapseRates.hh"

#include "pism/coupler/util/options.hh"
#include "pism/coupler/util/lapse_rates.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace atmosphere {

LapseRates::LapseRates(IceGrid::ConstPtr grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in),
  m_surface(grid, "ice_surface_elevation", WITHOUT_GHOSTS) {

  m_precip_lapse_rate = m_config->get_double("atmosphere.lapse_rate.precipitation_lapse_rate",
                                             "(m / s) / m");

  {
    options::Real T_lapse_rate("-temp_lapse_rate",
                               "Elevation lapse rate for the temperature, in K per km",
                               m_temp_lapse_rate);
    m_temp_lapse_rate = units::convert(m_sys, T_lapse_rate, "K/km", "K/m");
  }

  {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.lapse_rate");

    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_reference_surface = IceModelVec2T::ForcingField(m_grid,
                                                      file,
                                                      "usurf",
                                                      "", // no standard name
                                                      buffer_size,
                                                      evaluations_per_year,
                                                      periodic);
    m_reference_surface->set_attrs("climate_forcing", "ice surface elevation", "m",
                                   "surface_altitude", 0);
  }

  m_precipitation = allocate_precipitation(grid);
  m_temperature   = allocate_temperature(grid);
}

LapseRates::~LapseRates() {
  // empty
}

void LapseRates::init_impl(const Geometry &geometry) {
  using units::convert;

  m_input_model->init(geometry);

  m_log->message(2,
                 "   [using air temperature and precipitation lapse corrections]\n");

  m_log->message(2,
                 "   air temperature lapse rate: %3.3f K per km\n"
                 "   precipitation lapse rate:   %3.3f m year-1 per km\n",
                 convert(m_sys, m_temp_lapse_rate, "K / m", "K / km"),
                 convert(m_sys, m_precip_lapse_rate, "(m / s) / m", "(m / year) / km"));

  ForcingOptions opt(*m_grid->ctx(), "atmosphere.lapse_rate");

  m_reference_surface->init(opt.filename, opt.period, opt.reference_time);
}

void LapseRates::update_impl(const Geometry &geometry, double t, double dt) {

  m_input_model->update(geometry, t, dt);

  m_reference_surface->update(t, dt);
  m_reference_surface->interp(t + 0.5*dt);

  // make a copy of the surface elevation so that it is available in methods computing
  // temperature and precipitation time series
  m_surface.copy_from(geometry.ice_surface_elevation);

  // precipitation
  {
    m_precipitation->copy_from(m_input_model->mean_precipitation());

    lapse_rate_correction(m_surface, *m_reference_surface,
                          m_precip_lapse_rate, *m_precipitation);
  }

  // temperature
  {
    m_temperature->copy_from(m_input_model->mean_annual_temp());

    lapse_rate_correction(m_surface, *m_reference_surface,
                          m_temp_lapse_rate, *m_temperature);
  }
}

const IceModelVec2S& LapseRates::mean_annual_temp_impl() const {
  return *m_temperature;
}

const IceModelVec2S& LapseRates::mean_precipitation_impl() const {
  return *m_precipitation;
}

void LapseRates::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();

  m_reference_surface->begin_access();
  m_surface.begin_access();
}

void LapseRates::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();

  m_reference_surface->end_access();
  m_surface.end_access();
}

void LapseRates::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_reference_surface->init_interpolation(ts);
}

void LapseRates::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  std::vector<double> usurf(m_ts_times.size());

  m_input_model->temp_time_series(i, j, result);

  m_reference_surface->interp(i, j, usurf);

  for (unsigned int m = 0; m < m_ts_times.size(); ++m) {
    result[m] -= m_temp_lapse_rate * (m_surface(i, j) - usurf[m]);
  }
}

void LapseRates::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  std::vector<double> usurf(m_ts_times.size());

  m_input_model->precip_time_series(i, j, result);

  m_reference_surface->interp(i, j, usurf);

  for (unsigned int m = 0; m < m_ts_times.size(); ++m) {
    result[m] -= m_precip_lapse_rate * (m_surface(i, j) - usurf[m]);
  }
}

} // end of namespace atmosphere
} // end of namespace pism
