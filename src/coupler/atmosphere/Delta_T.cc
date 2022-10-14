// Copyright (C) 2011--2022 PISM Authors
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
#include "pism/util/ScalarForcing.hh"
#include "pism/coupler/util/options.hh"

namespace pism {
namespace atmosphere {

Delta_T::Delta_T(IceGrid::ConstPtr grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in) {

  std::string
    prefix         = "atmosphere.delta_T",
    variable_name  = "delta_T",
    long_name      = "near-surface air temperature offsets",
    units          = "Kelvin",
    external_units = "Kelvin";

  ForcingOptions opt(*m_grid->ctx(), prefix);

  // will be closed at the end of scope
  File input(m_grid->com, opt.filename, PISM_GUESS, PISM_READONLY);

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
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    m_2d_offsets = std::make_shared<array::Forcing>(m_grid,
                                                    input,
                                                    variable_name,
                                                    "", // no standard name
                                                    buffer_size,
                                                    opt.periodic);

    m_2d_offsets->set_attrs("climate_forcing",
                            long_name, units, external_units,
                            "", // no standard name
                            0);
  }

  m_temperature = allocate_temperature(grid);
}

void Delta_T::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing near-surface air temperature offsets...\n");

  if (m_2d_offsets) {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.delta_T");
    m_2d_offsets->init(opt.filename, opt.periodic);
  }
}

void Delta_T::init_timeseries_impl(const std::vector<double> &ts) const {
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

void Delta_T::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();

  if (m_2d_offsets) {
    m_2d_offsets->begin_access();
  }
}

void Delta_T::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();

  if (m_2d_offsets) {
    m_2d_offsets->end_access();
  }
}

void Delta_T::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);
  m_temperature->copy_from(m_input_model->air_temperature());

  if (m_1d_offsets) {
    m_temperature->shift(m_1d_offsets->value(t + 0.5 * dt));
  }

  if (m_2d_offsets) {
    m_2d_offsets->update(t, dt);
    m_2d_offsets->average(t, dt);

    auto &T = *m_temperature;
    const auto &delta = *m_2d_offsets;

    array::AccessScope list{&T, &delta};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      T(i, j) += delta(i, j);
    }
  }
}

const array::Scalar& Delta_T::air_temperature_impl() const {
  return *m_temperature;
}

void Delta_T::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->temp_time_series(i, j, result);

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
