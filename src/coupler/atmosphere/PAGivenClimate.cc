// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "PAGivenClimate.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"

namespace pism {

PAGivenClimate::PAGivenClimate(IceGrid &g, const Config &conf)
  : PGivenClimate<PAModifier,AtmosphereModel>(g, conf, NULL)
{
  option_prefix = "-atmosphere_given";
  air_temp      = NULL;
  precipitation = NULL;

  // Cannot call allocate_PAGivenClimate() here, because some surface
  // models do not use atmosphere models *and* this is the default
  // atmosphere model.

}

PAGivenClimate::~PAGivenClimate() {
  // empty
}

PetscErrorCode PAGivenClimate::allocate_PAGivenClimate() {
  PetscErrorCode ierr;

  // will be de-allocated by the parent's destructor
  precipitation = new IceModelVec2T;
  air_temp      = new IceModelVec2T;

  m_fields["precipitation"] = precipitation;
  m_fields["air_temp"]      = air_temp;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  ierr = air_temp->create(grid, "air_temp", false); CHKERRQ(ierr);
  ierr = precipitation->create(grid, "precipitation", false); CHKERRQ(ierr);

  ierr = air_temp->set_attrs("climate_forcing", "near-surface air temperature",
                             "Kelvin", ""); CHKERRQ(ierr);
  ierr = precipitation->set_attrs("climate_forcing", "ice-equivalent precipitation rate",
                                  "m s-1", ""); CHKERRQ(ierr);
  ierr = precipitation->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  precipitation->write_in_glaciological_units = true;

  return 0;
}

void PAGivenClimate::init(Vars &) {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, grid.com,
             "* Initializing the atmosphere model reading near-surface air temperature\n"
             "  and ice-equivalent precipitation from a file...\n");

  if (air_temp == NULL || precipitation == NULL) {
    allocate_PAGivenClimate();
  }

  air_temp->init(filename, bc_period, bc_reference_time);
  precipitation->init(filename, bc_period, bc_reference_time);

  // read time-independent data right away:
  if (air_temp->get_n_records() == 1 && precipitation->get_n_records() == 1) {
    update(grid.time->current(), 0); // dt is irrelevant
  }
}

void PAGivenClimate::update(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  // compute mean precipitation
  precipitation->average(m_t, m_dt);

  // Average so that the mean_annual_temp() may be reported correctly (at least
  // in the "-surface pdd" case).
  air_temp->average(m_t, m_dt);
}

void PAGivenClimate::mean_precipitation(IceModelVec2S &result) {
  precipitation->copy_to(result);
}

void PAGivenClimate::mean_annual_temp(IceModelVec2S &result) {
  air_temp->copy_to(result);
}

void PAGivenClimate::temp_snapshot(IceModelVec2S &result) {
  air_temp->copy_to(result);
}

void PAGivenClimate::begin_pointwise_access() {

  air_temp->begin_access();
  precipitation->begin_access();
}

void PAGivenClimate::end_pointwise_access() {

  air_temp->end_access();
  precipitation->end_access();
}

void PAGivenClimate::temp_time_series(int i, int j, std::vector<double> &result) {

  air_temp->interp(i, j, result);
}

void PAGivenClimate::precip_time_series(int i, int j, std::vector<double> &result) {

  precipitation->interp(i, j, result);
}

void PAGivenClimate::init_timeseries(const std::vector<double> &ts) {

  air_temp->init_interpolation(ts);

  precipitation->init_interpolation(ts);

  m_ts_times = ts;
}


} // end of namespace pism
