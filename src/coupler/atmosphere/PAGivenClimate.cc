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
  air_temp      = NULL;
  precipitation = NULL;

  // will be de-allocated by the parent's destructor
  precipitation = new IceModelVec2T;
  air_temp      = new IceModelVec2T;

  m_fields["precipitation"] = precipitation;
  m_fields["air_temp"]      = air_temp;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  air_temp->create(m_grid, "air_temp");
  precipitation->create(m_grid, "precipitation");

  air_temp->set_attrs("climate_forcing", "near-surface air temperature",
                      "Kelvin", "");
  air_temp->metadata().set_double("valid_min", 0.0);
  air_temp->metadata().set_double("valid_max", 323.15); // 50 C

  precipitation->set_attrs("climate_forcing", "ice-equivalent precipitation rate",
                           "kg m-2 second-1", "");
  precipitation->metadata().set_string("glaciological_units", "kg m-2 year-1");
  precipitation->write_in_glaciological_units = true;
}

Given::~Given() {
  // empty
}

void Given::init() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the atmosphere model reading near-surface air temperature\n"
             "  and ice-equivalent precipitation from a file...\n");

  air_temp->init(m_filename, m_bc_period, m_bc_reference_time);
  precipitation->init(m_filename, m_bc_period, m_bc_reference_time);

  // read time-independent data right away:
  if (air_temp->get_n_records() == 1 && precipitation->get_n_records() == 1) {
    update(m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  // compute mean precipitation
  precipitation->average(m_t, m_dt);

  // Average so that the mean_annual_temp() may be reported correctly (at least
  // in the "-surface pdd" case).
  air_temp->average(m_t, m_dt);
}

void Given::mean_precipitation(IceModelVec2S &result) {
  result.copy_from(*precipitation);
}

void Given::mean_annual_temp(IceModelVec2S &result) {
  result.copy_from(*air_temp);
}

void Given::temp_snapshot(IceModelVec2S &result) {
  result.copy_from(*air_temp);
}

void Given::begin_pointwise_access() {

  air_temp->begin_access();
  precipitation->begin_access();
}

void Given::end_pointwise_access() {

  air_temp->end_access();
  precipitation->end_access();
}

void Given::temp_time_series(int i, int j, std::vector<double> &result) {

  air_temp->interp(i, j, result);
}

void Given::precip_time_series(int i, int j, std::vector<double> &result) {

  precipitation->interp(i, j, result);
}

void Given::init_timeseries(const std::vector<double> &ts) {

  air_temp->init_interpolation(ts);

  precipitation->init_interpolation(ts);

  m_ts_times = ts;
}


} // end of namespace atmosphere
} // end of namespace pism
