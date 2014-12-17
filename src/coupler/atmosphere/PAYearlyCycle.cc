// Copyright (C) 2008-2014 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

#include "PAYearlyCycle.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"
#include <stdexcept>

namespace pism {

PAYearlyCycle::PAYearlyCycle(IceGrid &g)
  : AtmosphereModel(g),
    m_air_temp_snapshot(g.config.get_unit_system(), "air_temp_snapshot", g) {

  m_snow_temp_july_day = m_config.get("snow_temp_july_day");

  // Allocate internal IceModelVecs:
  m_air_temp_mean_annual.create(m_grid, "air_temp_mean_annual", WITHOUT_GHOSTS);
  m_air_temp_mean_annual.set_attrs("diagnostic",
                                   "mean annual near-surface air temperature (without sub-year time-dependence or forcing)",
                                   "K", 
                                   "");  // no CF standard_name ??
  m_air_temp_mean_annual.metadata().set_string("source", m_reference);

  m_air_temp_mean_july.create(m_grid, "air_temp_mean_july", WITHOUT_GHOSTS);
  m_air_temp_mean_july.set_attrs("diagnostic",
                                 "mean July near-surface air temperature (without sub-year time-dependence or forcing)",
                                 "Kelvin",
                                 "");  // no CF standard_name ??
  m_air_temp_mean_july.metadata().set_string("source", m_reference);

  m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
  m_precipitation.set_attrs("climate_state", 
                            "mean annual ice-equivalent precipitation rate",
                            "m s-1", 
                            ""); // no CF standard_name ??
  m_precipitation.set_glaciological_units("m year-1");
  m_precipitation.write_in_glaciological_units = true;
  m_precipitation.set_time_independent(true);

  m_air_temp_snapshot.set_string("pism_intent", "diagnostic");
  m_air_temp_snapshot.set_string("long_name",
                         "snapshot of the near-surface air temperature");
  m_air_temp_snapshot.set_units("K");
}

PAYearlyCycle::~PAYearlyCycle() {
  // empty
}

//! Allocates memory and reads in the precipitaion data.
void PAYearlyCycle::init() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  bool do_regrid = false;
  int start = -1;
  find_pism_input(m_precip_filename, do_regrid, start);

  init_internal(m_precip_filename, do_regrid, start);
}

//! Read precipitation data from a given file.
void PAYearlyCycle::init_internal(const std::string &input_filename, bool do_regrid,
                                            unsigned int start) {
  // read precipitation rate from file
  verbPrintf(2, m_grid.com,
             "    reading mean annual ice-equivalent precipitation rate 'precipitation'\n"
             "      from %s ... \n",
             input_filename.c_str());
  if (do_regrid == true) {
    m_precipitation.regrid(input_filename, CRITICAL); // fails if not found!
  } else {
    m_precipitation.read(input_filename, start); // fails if not found!
  }
}

void PAYearlyCycle::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  result.insert("precipitation");

  if (keyword == "big") {
    result.insert("air_temp_mean_annual");
    result.insert("air_temp_mean_july");
    result.insert("air_temp_snapshot");
  }
}


void PAYearlyCycle::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  if (set_contains(vars, "air_temp_snapshot")) {
    m_air_temp_snapshot.define(nc, nctype, false);
  }

  if (set_contains(vars, "air_temp_mean_annual")) {
    m_air_temp_mean_annual.define(nc, nctype);
  }

  if (set_contains(vars, "air_temp_mean_july")) {
    m_air_temp_mean_july.define(nc, nctype);
  }

  if (set_contains(vars, "precipitation")) {
    m_precipitation.define(nc, nctype);
  }
}


void PAYearlyCycle::write_variables(const std::set<std::string> &vars, const PIO &nc) {

  if (set_contains(vars, "air_temp_snapshot")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "air_temp_snapshot", WITHOUT_GHOSTS);
    tmp.metadata() = m_air_temp_snapshot;

    temp_snapshot(tmp);

    tmp.write(nc);
  }

  if (set_contains(vars, "air_temp_mean_annual")) {
    m_air_temp_mean_annual.write(nc);
  }

  if (set_contains(vars, "air_temp_mean_july")) {
    m_air_temp_mean_july.write(nc);
  }

  if (set_contains(vars, "precipitation")) {
    m_precipitation.write(nc);
  }
}

//! Copies the stored precipitation field into result.
void PAYearlyCycle::mean_precipitation(IceModelVec2S &result) {
  m_precipitation.copy_to(result);
}

//! Copies the stored mean annual near-surface air temperature field into result.
void PAYearlyCycle::mean_annual_temp(IceModelVec2S &result) {
  m_air_temp_mean_annual.copy_to(result);
}

void PAYearlyCycle::init_timeseries(const std::vector<double> &ts) {
  // constants related to the standard yearly cycle
  const double
    julyday_fraction = m_grid.time->day_of_the_year_to_day_fraction(m_snow_temp_july_day);

  size_t N = ts.size();

  m_ts_times.resize(N);
  m_cosine_cycle.resize(N);
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    double tk = m_grid.time->year_fraction(ts[k]) - julyday_fraction;

    m_ts_times[k] = ts[k];
    m_cosine_cycle[k] = cos(2.0 * M_PI * tk);
  }
}

void PAYearlyCycle::precip_time_series(int i, int j, std::vector<double> &result) {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

void PAYearlyCycle::temp_time_series(int i, int j, std::vector<double> &result) {

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] = m_air_temp_mean_annual(i,j) + (m_air_temp_mean_july(i,j) - m_air_temp_mean_annual(i,j)) * m_cosine_cycle[k];
  }
}

void PAYearlyCycle::temp_snapshot(IceModelVec2S &result) {
  const double
    julyday_fraction = m_grid.time->day_of_the_year_to_day_fraction(m_snow_temp_july_day),
    T                = m_grid.time->year_fraction(m_t + 0.5 * m_dt) - julyday_fraction,
    cos_T            = cos(2.0 * M_PI * T);

  IceModelVec::AccessList list;
  list.add(result);
  list.add(m_air_temp_mean_annual);
  list.add(m_air_temp_mean_july);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    result(i,j) = m_air_temp_mean_annual(i,j) + (m_air_temp_mean_july(i,j) - m_air_temp_mean_annual(i,j)) * cos_T;
  }
}

void PAYearlyCycle::begin_pointwise_access() {
  m_air_temp_mean_annual.begin_access();
  m_air_temp_mean_july.begin_access();
  m_precipitation.begin_access();
}

void PAYearlyCycle::end_pointwise_access() {
  m_air_temp_mean_annual.end_access();
  m_air_temp_mean_july.end_access();
  m_precipitation.end_access();
}


} // end of namespace pism
