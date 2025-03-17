// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023, 2024, 2025 PISM Authors
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

#include "pism/coupler/atmosphere/GivenClimate.hh"

#include "pism/coupler/util/options.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Config.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace atmosphere {

Given::Given(std::shared_ptr<const Grid> g)
  : AtmosphereModel(g, std::shared_ptr<AtmosphereModel>()) {
  ForcingOptions opt(*m_grid->ctx(), "atmosphere.given");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    auto interp_type = m_config->get_string("atmosphere.given.air_temperature_interpolation");

    InterpolationType interpolation = interp_type == "piecewise_linear" ? LINEAR : PIECEWISE_CONSTANT;

    m_air_temp = std::make_shared<array::Forcing>(m_grid,
                                             file,
                                             "air_temp",
                                             "", // no standard name
                                             buffer_size,
                                             opt.periodic,
                                             interpolation);

    m_precipitation = std::make_shared<array::Forcing>(m_grid,
                                                       file,
                                                       "precipitation",
                                                       "", // no standard name
                                                       buffer_size,
                                                       opt.periodic);
  }

  {
    m_air_temp->metadata(0)
        .long_name("mean annual near-surface air temperature")
        .units("kelvin");
    m_air_temp->metadata(0)["valid_range"] = { 0.0, 323.15 }; // (0 C, 50 C)
  }
  {
    m_precipitation->metadata(0)
        .long_name("precipitation rate")
        .units("kg m^-2 second^-1")
        .output_units("kg m^-2 year^-1")
        .standard_name("precipitation_flux");
  }
}

void Given::init_impl(const Geometry &geometry) {
  m_log->message(2,
             "* Initializing the atmosphere model reading near-surface air temperature\n"
             "  and ice-equivalent precipitation from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "atmosphere.given");

  m_air_temp->init(opt.filename, opt.periodic);
  m_precipitation->init(opt.filename, opt.periodic);

  // read time-independent data right away:
  if (m_air_temp->buffer_size() == 1 && m_precipitation->buffer_size() == 1) {
    update(geometry, time().current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;

  m_precipitation->update(t, dt);
  m_air_temp->update(t, dt);

  m_precipitation->average(t, dt);
  m_air_temp->average(t, dt);
}

const array::Scalar& Given::precipitation_impl() const {
  return *m_precipitation;
}

const array::Scalar& Given::air_temperature_impl() const {
  return *m_air_temp;
}

void Given::begin_pointwise_access_impl() const {

  m_air_temp->begin_access();
  m_precipitation->begin_access();
}

void Given::end_pointwise_access_impl() const {

  m_air_temp->end_access();
  m_precipitation->end_access();
}

void Given::temp_time_series_impl(int i, int j, std::vector<double> &result) const {

  m_air_temp->interp(i, j, result);
}

void Given::precip_time_series_impl(int i, int j, std::vector<double> &result) const {

  m_precipitation->interp(i, j, result);
}

void Given::init_timeseries_impl(const std::vector<double> &ts) const {

  m_air_temp->init_interpolation(ts);

  m_precipitation->init_interpolation(ts);

  m_ts_times = ts;
}


} // end of namespace atmosphere
} // end of namespace pism
