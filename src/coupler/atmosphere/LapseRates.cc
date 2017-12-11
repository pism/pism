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

#include <cassert>
#include <gsl/gsl_math.h>

#include "LapseRates.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace atmosphere {

LapseRates::LapseRates(IceGrid::ConstPtr g, AtmosphereModel* in)
  : PLapseRates<AtmosphereModel,PAModifier>(g, in) {
  m_precip_lapse_rate = 0.0;
  m_option_prefix     = "-atmosphere_lapse_rate";
  m_surface = NULL;

  m_precip_scale_factor = 0.0;
  do_precip_scale  = m_config->get_boolean("atmosphere.precip_lapse_scaling");

}

LapseRates::~LapseRates() {
  // empty
}

void LapseRates::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2,
             "  [using air temperature and precipitation lapse corrections]\n");

  init_internal();

  m_precip_lapse_rate  = m_config->get_double("atmosphere.precip_lapse_rate");

  //This is basically temperature lapse rate 8.2 K/Km as in TemperaturPIK times precipitation scale rate 5%/K )
  m_precip_scale_factor  = m_config->get_double("atmosphere.precip_lapse_scale_factor");

  if (do_precip_scale){

    m_log->message(2,
              "   air temperature lapse rate: %3.3f K per km\n"
              "   precipitation scale factor: %3.3f per km\n",
              m_temp_lapse_rate, m_precip_scale_factor);
  } else {

    m_log->message(2,
             "   air temperature lapse rate: %3.3f K per km\n"
             "   precipitation lapse rate:   %3.3f m year-1 per km\n",
             m_temp_lapse_rate, m_precip_lapse_rate);
  }

  m_temp_lapse_rate = units::convert(m_sys, m_temp_lapse_rate, "K/km", "K/m");

  m_precip_lapse_rate = units::convert(m_sys, m_precip_lapse_rate, "m year-1 / km", "m second-1 / m");

  m_precip_scale_factor = units::convert(m_sys, m_precip_scale_factor, "km-1", "m-1");

  m_surface = m_grid->variables().get_2d_scalar("surface_altitude");
}

void LapseRates::mean_precipitation_impl(IceModelVec2S &result) const {
  m_input_model->mean_precipitation(result);

  if (do_precip_scale) {
    lapse_rate_scale(result, m_precip_scale_factor);
  } else {
    lapse_rate_correction(result, m_precip_lapse_rate);
  }
}

void LapseRates::mean_annual_temp_impl(IceModelVec2S &result) const {
  m_input_model->mean_annual_temp(result);
  lapse_rate_correction(result, m_temp_lapse_rate);
}

void LapseRates::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();
  m_reference_surface.begin_access();

  assert(m_surface != NULL);
  m_surface->begin_access();
}

void LapseRates::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();
  m_reference_surface.end_access();

  assert(m_surface != NULL);
  m_surface->end_access();
}

void LapseRates::init_timeseries_impl(const std::vector<double> &ts) const {
  PAModifier::init_timeseries_impl(ts);

  m_reference_surface.init_interpolation(ts);

}

void LapseRates::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  std::vector<double> usurf(m_ts_times.size());

  m_input_model->temp_time_series(i, j, result);

  m_reference_surface.interp(i, j, usurf);

  assert(m_surface != NULL);

  for (unsigned int m = 0; m < m_ts_times.size(); ++m) {
    result[m] -= m_temp_lapse_rate * ((*m_surface)(i, j) - usurf[m]);
  }
}

void LapseRates::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  std::vector<double> usurf(m_ts_times.size());

  m_input_model->precip_time_series(i, j, result);

  m_reference_surface.interp(i, j, usurf);

  assert(m_surface != NULL);

  for (unsigned int m = 0; m < m_ts_times.size(); ++m) {

    if (do_precip_scale) {
      result[m] *= exp(-m_precip_scale_factor * ((*m_surface)(i, j) - usurf[m]));
    } else {
      result[m] -= m_precip_lapse_rate * ((*m_surface)(i, j) - usurf[m]);
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
