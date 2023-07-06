// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023 PISM Authors
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

#include "pism/coupler/atmosphere/ElevationChange.hh"

#include <cmath>                // std::exp()

#include "pism/coupler/util/options.hh"
#include "pism/coupler/util/lapse_rates.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace atmosphere {

ElevationChange::ElevationChange(std::shared_ptr<const Grid> grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in),
  m_surface(grid, "ice_surface_elevation") {

  m_precip_lapse_rate = m_config->get_number("atmosphere.elevation_change.precipitation.lapse_rate",
                                             "(kg m-2 / s) / m");

  m_precip_temp_lapse_rate = m_config->get_number("atmosphere.elevation_change.precipitation.temp_lapse_rate",
                                                  "K / m");
  m_precip_exp_factor = m_config->get_number("atmosphere.precip_exponential_factor_for_temperature");

  m_temp_lapse_rate = m_config->get_number("atmosphere.elevation_change.temperature_lapse_rate",
                                           "K / m");

  {
    auto method = m_config->get_string("atmosphere.elevation_change.precipitation.method");
    m_precip_method = method == "scale" ? SCALE : SHIFT;
  }

  {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.elevation_change");

    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_reference_surface = std::make_shared<array::Forcing>(m_grid,
                                                      file,
                                                      "usurf",
                                                      "", // no standard name
                                                      buffer_size,
                                                      opt.periodic,
                                                      LINEAR);
    m_reference_surface->metadata()
        .intent("climate_forcing")
        .long_name("ice surface elevation")
        .units("m")
        .standard_name("surface_altitude");
  }

  m_precipitation = allocate_precipitation(grid);
  m_temperature   = allocate_temperature(grid);
}

void ElevationChange::init_impl(const Geometry &geometry) {
  using units::convert;

  m_input_model->init(geometry);

  m_log->message(2,
                 "   [using elevation-change-dependent adjustments of air temperature and precipitation]\n");

    m_log->message(2,
                   "   air temperature lapse rate: %3.3f K per km\n",
                   convert(m_sys, m_temp_lapse_rate, "K / m", "K / km"));

  if (m_precip_method == SHIFT) {
    m_log->message(2,
                   "   precipitation lapse rate:   %3.3f (kg m-2 year-1) per km\n",
                   convert(m_sys, m_precip_lapse_rate, "(kg m-2 / s) / m", "(kg m-2 / year) / km"));
  } else {
    m_log->message(2,
                   "   precipitation scaling factor with temperature: %3.3f Kelvin-1\n"
                   "   temperature lapse rate: %3.3f K per km\n",
                   m_precip_exp_factor,
                   convert(m_sys, m_precip_temp_lapse_rate, "K / m", "K / km"));
  }

  ForcingOptions opt(*m_grid->ctx(), "atmosphere.elevation_change");

  m_reference_surface->init(opt.filename, opt.periodic);
}

void ElevationChange::update_impl(const Geometry &geometry, double t, double dt) {

  m_input_model->update(geometry, t, dt);

  m_reference_surface->update(t, dt);
  m_reference_surface->interp(t + 0.5*dt);

  // make a copy of the surface elevation so that it is available in methods computing
  // temperature and precipitation time series
  m_surface.copy_from(geometry.ice_surface_elevation);

  const auto &reference_surface = *m_reference_surface;

  // temperature
  {
    m_temperature->copy_from(m_input_model->air_temperature());

    lapse_rate_correction(m_surface, reference_surface,
                          m_temp_lapse_rate, *m_temperature);
  }

  // precipitation
  {
    m_precipitation->copy_from(m_input_model->precipitation());

    switch (m_precip_method) {
    case SCALE:
      {
        array::AccessScope list{&m_surface, &reference_surface, m_precipitation.get()};

        for (auto p = m_grid->points(); p; p.next()) {
          const int i = p.i(), j = p.j();

          double dT = -m_precip_temp_lapse_rate * (m_surface(i, j) - reference_surface(i, j));

          (*m_precipitation)(i, j) *= std::exp(m_precip_exp_factor * dT);
        }

      }
      break;
    case SHIFT:
    default:
      {
        lapse_rate_correction(m_surface, *m_reference_surface,
                              m_precip_lapse_rate, *m_precipitation);
      }
      break;
    }
  }
}

const array::Scalar& ElevationChange::air_temperature_impl() const {
  return *m_temperature;
}

const array::Scalar& ElevationChange::precipitation_impl() const {
  return *m_precipitation;
}

void ElevationChange::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();

  m_reference_surface->begin_access();
  m_surface.begin_access();
}

void ElevationChange::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();

  m_reference_surface->end_access();
  m_surface.end_access();
}

void ElevationChange::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_reference_surface->init_interpolation(ts);
}

void ElevationChange::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  std::vector<double> reference_surface(m_ts_times.size());

  m_input_model->temp_time_series(i, j, result);

  m_reference_surface->interp(i, j, reference_surface);

  for (unsigned int m = 0; m < m_ts_times.size(); ++m) {
    result[m] -= m_temp_lapse_rate * (m_surface(i, j) - reference_surface[m]);
  }
}

void ElevationChange::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  auto N = m_ts_times.size();
  std::vector<double> reference_surface(N);

  m_input_model->precip_time_series(i, j, result);

  m_reference_surface->interp(i, j, reference_surface);

  switch (m_precip_method) {
  case SCALE:
    {
      for (unsigned int m = 0; m < N; ++m) {
        double dT = -m_precip_temp_lapse_rate * (m_surface(i, j) - reference_surface[m]);
        result[m] *= std::exp(m_precip_exp_factor * dT);
      }
    }
    break;
  case SHIFT:
    for (unsigned int m = 0; m < N; ++m) {
      result[m] -= m_precip_lapse_rate * (m_surface(i, j) - reference_surface[m]);
    }
    break;
  }
}

} // end of namespace atmosphere
} // end of namespace pism
