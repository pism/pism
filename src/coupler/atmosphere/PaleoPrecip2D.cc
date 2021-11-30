// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2021 PISM Authors
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

#include "PaleoPrecip2D.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/io_helpers.hh"

namespace pism {
namespace atmosphere {

PaleoPrecip2D::PaleoPrecip2D(IceGrid::ConstPtr g, AtmosphereModel* in)
  : PGivenClimate<PAModifier,AtmosphereModel>(g, in) {
  m_option_prefix  = "-atmosphere_paleo_precip_2d";

  process_options();

  m_precipexpfactor = m_config->get_double("atmosphere.precip_exponential_factor_for_temperature");

  // will be de-allocated by the parent's destructor
  {
    m_delta_T = new IceModelVec2T;
    m_fields["delta_T"] = m_delta_T;
  }

  // calls set_n_records(), so ...->create() has to be called after set_vec_parameters()
  set_vec_parameters({});

  {
    m_delta_T->create(m_grid, "delta_T");
    m_delta_T->set_attrs("climate_forcing",
                         "anomaly of the near-surface air temperature",
                         "Kelvin", "");
  }
}

PaleoPrecip2D::~PaleoPrecip2D()
{
  // empty
}

void PaleoPrecip2D::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2,
             "* Initializing the -atmosphere ...,anomaly code...\n");

  m_log->message(2,
             "    reading anomalies from %s ...\n",
             m_filename.c_str());

  m_delta_T->init(m_filename, m_bc_period, m_bc_reference_time);
}

void PaleoPrecip2D::update_impl(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  m_delta_T->average(m_t, m_dt);
}

void PaleoPrecip2D::mean_precipitation_impl(IceModelVec2S &result) const {
  m_input_model->mean_precipitation(result);

  IceModelVec::AccessList list{&result, m_delta_T};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double dT = (*m_delta_T)(i, j);

    result(i, j) *= exp(m_precipexpfactor * dT);
  }
}

void PaleoPrecip2D::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();
  m_delta_T->begin_access();
}

void PaleoPrecip2D::end_pointwise_access_impl() const {
  m_input_model->end_pointwise_access();
  m_delta_T->end_access();
}

void PaleoPrecip2D::init_timeseries_impl(const std::vector<double> &ts) const {
  PAModifier::init_timeseries_impl(ts);

  m_delta_T->init_interpolation(ts);
}

void PaleoPrecip2D::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->precip_time_series(i, j, result);

  m_temp_anomaly.reserve(m_ts_times.size());
  m_delta_T->interp(i, j, m_temp_anomaly);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] *= exp(m_precipexpfactor * m_temp_anomaly[k]);
  }
}

} // end of namespace atmosphere
} // end of namespace pism
