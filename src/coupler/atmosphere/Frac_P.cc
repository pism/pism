// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2021 PISM Authors
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

#include "Frac_P.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/util/io/File.hh"
#include "pism/coupler/util/options.hh"

namespace pism {
namespace atmosphere {

Frac_P::Frac_P(IceGrid::ConstPtr grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in) {

  std::string
    prefix        = "atmosphere.frac_P",
    variable_name = "frac_P",
    long_name     = "precipitation multiplier, pure fraction",
    units         = "1";

  ForcingOptions opt(*m_grid->ctx(), prefix);

  // will be closed at the end of scope
  File input(m_grid->com, opt.filename, PISM_NETCDF3, PISM_READONLY);

  // Assume that we are expected to use 1D scaling if the input file contains a scalar
  // time-series.
  bool scalar = input.dimensions(variable_name).size() == 1;

  if (scalar) {
    m_1d_scaling.reset(new ScalarForcing(*grid->ctx(),
                                         prefix,
                                         variable_name,
                                         units, units,
                                         long_name));
  } else {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_number("input.forcing.evaluations_per_year");

    m_2d_scaling = IceModelVec2T::ForcingField(m_grid,
                                               input,
                                               variable_name,
                                               "", // no standard name
                                               buffer_size,
                                               evaluations_per_year,
                                               opt.periodic);

    m_2d_scaling->set_attrs("climate_forcing",
                            long_name, units, units, "", 0);
  }

  m_precipitation = allocate_precipitation(grid);
}

Frac_P::~Frac_P() {
  // empty
}

void Frac_P::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2, "* Initializing precipitation forcing using scalar multipliers...\n");

  if (m_2d_scaling) {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.frac_P");
    m_2d_scaling->init(opt.filename, opt.periodic);
  }
}

void Frac_P::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);

  m_scaling_values.resize(ts.size());

  if (m_1d_scaling) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      m_scaling_values[k] = m_1d_scaling->value(ts[k]);
    }
  }

  if (m_2d_scaling) {
    m_2d_scaling->init_interpolation(ts);
  }
}

void Frac_P::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();

  if (m_2d_scaling) {
    m_2d_scaling->begin_access();
  }
}

void Frac_P::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();

  if (m_2d_scaling) {
    m_2d_scaling->end_access();
  }
}

void Frac_P::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);
  m_precipitation->copy_from(m_input_model->mean_precipitation());

  if (m_1d_scaling) {
    m_precipitation->scale(m_1d_scaling->value(t + 0.5 * dt));
  }

  if (m_2d_scaling) {
    m_2d_scaling->update(t, dt);
    m_2d_scaling->average(t, dt);

    IceModelVec2S &P = *m_precipitation;
    IceModelVec2T &S = *m_2d_scaling;

    IceModelVec::AccessList list{&P, &S};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      P(i, j) *= S(i, j);
    }
  }
}

const IceModelVec2S& Frac_P::mean_precipitation_impl() const {
  return *m_precipitation;
}

void Frac_P::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->precip_time_series(i, j, result);

  if (m_2d_scaling) {
    // m_scaling_values was resized in init_interpolation and so it should have enough
    // elements
    m_2d_scaling->interp(i, j, m_scaling_values);
  } else if (m_1d_scaling) {
    // empty: m_scaling_values were set in init_timeseries_impl()
  }

  for (size_t k = 0; k < result.size(); ++k) {
    result[k] *= m_scaling_values[k];
  }
}

} // end of namespace atmosphere
} // end of namespace pism
