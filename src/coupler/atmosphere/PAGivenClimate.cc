// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include "PAGivenClimate.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"

namespace pism {
namespace atmosphere {

Given::Given(IceGrid::ConstPtr g)
  : PGivenClimate<PAModifier,AtmosphereModel>(g, NULL)
{
  m_option_prefix = "-atmosphere_given";
  m_air_temp      = NULL;
  m_precipitation = NULL;

  // will be de-allocated by the parent's destructor
  m_precipitation = new IceModelVec2T;
  m_air_temp      = new IceModelVec2T;

  m_fields["precipitation"] = m_precipitation;
  m_fields["air_temp"]      = m_air_temp;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  {
    m_air_temp->create(m_grid, "air_temp");
    m_air_temp->set_attrs("diagnostic", "mean annual near-surface air temperature",
                          "Kelvin", "", 0);
    m_air_temp->metadata(0).set_double("valid_min", 0.0);
    m_air_temp->metadata(0).set_double("valid_max", 323.15); // 50 C
  }
  {
    m_precipitation->create(m_grid, "precipitation");
    m_precipitation->set_attrs("model_state", "precipitation rate",
                               "kg m-2 second-1", "", 0);
    m_precipitation->metadata(0).set_string("glaciological_units", "kg m-2 year-1");
    m_precipitation->write_in_glaciological_units = true;
  }
}

Given::~Given() {
  // empty
}

void Given::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the atmosphere model reading near-surface air temperature\n"
             "  and ice-equivalent precipitation from a file...\n");

  m_air_temp->init(m_filename, m_bc_period, m_bc_reference_time);
  m_precipitation->init(m_filename, m_bc_period, m_bc_reference_time);

  // read time-independent data right away:
  if (m_air_temp->get_n_records() == 1 && m_precipitation->get_n_records() == 1) {
    update(m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  // compute mean precipitation
  m_precipitation->average(m_t, m_dt);

  // Average so that the mean_annual_temp() may be reported correctly (at least
  // in the "-surface pdd" case).
  m_air_temp->average(m_t, m_dt);
}

void Given::mean_precipitation_impl(IceModelVec2S &result) const {
  result.copy_from(*m_precipitation);
}

void Given::mean_annual_temp_impl(IceModelVec2S &result) const {
  result.copy_from(*m_air_temp);
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
