// Copyright (C) 2011--2024 PISM Authors
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

#include "pism/coupler/atmosphere/Delta_P.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace atmosphere {

Delta_P::Delta_P(std::shared_ptr<const Grid> grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, std::move(in)) {

  std::string
    prefix         = "atmosphere.delta_P",
    variable_name  = "delta_P",
    long_name      = "precipitation offsets",
    units          = "kg m^-2 second^-1",
    external_units = "kg m^-2 year^-1";

  ForcingOptions opt(*m_grid->ctx(), prefix);

  // will be closed at the end of scope
  File input(m_grid->com, opt.filename, io::PISM_GUESS, io::PISM_READONLY);

  // Assume that we are expected to use 1D scaling if the input file contains a scalar
  // time-series.
  bool scalar = input.dimensions(variable_name).size() == 1;

  if (scalar) {
    m_1d_offsets.reset(new ScalarForcing(*grid->ctx(),
                                         prefix,
                                         variable_name,
                                         units, external_units,
                                         long_name));
  } else {
    auto buffer_size = m_config->get_number("input.forcing.buffer_size");

    m_2d_offsets = std::make_shared<array::Forcing>(m_grid,
                                                    input,
                                                    variable_name,
                                                    "", // no standard name
                                                    static_cast<unsigned int>(buffer_size),
                                                    opt.periodic);

    m_2d_offsets->metadata()
        .long_name(long_name)
        .units(units)
        .output_units(external_units);
  }

  m_precipitation = allocate_precipitation(grid);
}

void Delta_P::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing precipitation forcing using scalar offsets...\n");

  if (m_2d_offsets) {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.delta_P");
    m_2d_offsets->init(opt.filename, opt.periodic);
  }

}

void Delta_P::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_offset_values.resize(ts.size());

  if (m_1d_offsets) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      m_offset_values[k] = m_1d_offsets->value(ts[k]);
    }
  }

  if (m_2d_offsets) {
    m_2d_offsets->init_interpolation(ts);
  }
}

void Delta_P::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();

  if (m_2d_offsets) {
    m_2d_offsets->begin_access();
  }
}

void Delta_P::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();

  if (m_2d_offsets) {
    m_2d_offsets->end_access();
  }
}

void Delta_P::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);
  m_precipitation->copy_from(m_input_model->precipitation());

  if (m_1d_offsets) {
    m_precipitation->shift(m_1d_offsets->value(t + 0.5 * dt));
  }

  if (m_2d_offsets) {
    m_2d_offsets->update(t, dt);
    m_2d_offsets->average(t, dt);

    m_precipitation->add(1.0, *m_2d_offsets);
  }
}

const array::Scalar& Delta_P::precipitation_impl() const {
  return *m_precipitation;
}

void Delta_P::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->precip_time_series(i, j, result);

  if (m_2d_offsets) {
    // m_offset_values was resized in init_interpolation and so it should have enough
    // elements
    m_2d_offsets->interp(i, j, m_offset_values);
  } else if (m_1d_offsets) {
    // empty: m_offset_values were set in init_timeseries_impl()
  }

  for (size_t k = 0; k < result.size(); ++k) {
    result[k] += m_offset_values[k];
  }
}

} // end of namespace atmosphere
} // end of namespace pism
