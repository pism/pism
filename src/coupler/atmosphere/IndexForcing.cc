// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#include "pism/util/pism_options.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/Time.hh"
#include "pism/geometry/Geometry.hh"

#include "IndexForcing.hh"

namespace pism {
namespace atmosphere {

IndexForcing::IndexForcing(IceGrid::ConstPtr g)
  : AtmosphereModel(g, nullptr) {

  m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
  m_precipitation.set_attrs("model_state", "precipitation rate",
                            "kg m-2 second-1",
                            "precipitation_flux", 0);
  m_precipitation.metadata(0).set_string("glaciological_units", "kg m-2 year-1");

  m_air_temp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
  m_air_temp.set_attrs("model_state", "mean annual near-surface air temperature",
                           "Kelvin", "", 0);

  m_h0.create(m_grid, "usurf_0", WITHOUT_GHOSTS);
  m_h0.set_attrs("climate_state", "surface elevation at t0", "m", "");
  m_h0.set_time_independent(true);

  m_h1.create(m_grid, "usurf_1", WITHOUT_GHOSTS);
  m_h1.set_attrs("climate_state", "surface elevation at t1", "m", "");
  m_h1.set_time_independent(true);

  {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.index");

    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_index = new Timeseries(*m_grid, "glac_index", m_config->get_string("time.dimension_name"));
    m_index->variable().set_string("units", "1");
    m_index->variable().set_string("long_name", "glacial index");
    m_index->dimension().set_string("units", m_grid->ctx()->time()->units_string());
    m_index_reference_time = opt.reference_time;
    m_index_period = opt.period;

    m_T0 = IceModelVec2T::ForcingField(m_grid,
                                      file,
                                      "airtemp_0",
                                      "", // no standard name
                                      buffer_size,
                                      evaluations_per_year,
                                      true);
    m_T0->set_attrs("climate_forcing", "air temperature at t0", "K", "");

    m_T1 = IceModelVec2T::ForcingField(m_grid,
                                      file,
                                      "airtemp_1",
                                      "", // no standard name
                                      buffer_size,
                                      evaluations_per_year,
                                      true);
    m_T1->set_attrs("climate_forcing", "air temperature at t1", "K", "");

    m_P0 = IceModelVec2T::ForcingField(m_grid,
                                      file,
                                      "precip_0",
                                      "", // no standard name
                                      buffer_size,
                                      evaluations_per_year,
                                      true);
    m_P0->set_attrs("climate_forcing", "precipitation at t0", "kg m-2 second-1", "");

    m_P1 = IceModelVec2T::ForcingField(m_grid,
                                      file,
                                      "precip_1",
                                      "", // no standard name
                                      buffer_size,
                                      evaluations_per_year,
                                      true);
    m_P1->set_attrs("climate_forcing", "precipitation at t1", "kg m-2 second-1", "");
  }

  m_temp_lapse_rate = 5.; // temp lapse rate [K/km]
  m_precip_decay_rate = log(2) / 1.; // precip_decay rate [1/km]
  m_precip_thresh_height = 5. ; // threshold height for precip [km]
}

IndexForcing::~IndexForcing() {
  // Do nothing....
}

void IndexForcing::init_impl(const Geometry &geometry) {
  m_temp_lapse_rate = options::Real("-temp_lapse_rate",
                                    "Elevation lapse rate for the temperature,"
                                    "in K per km",
                                    m_temp_lapse_rate);
  m_temp_lapse_rate = units::convert(m_sys, m_temp_lapse_rate, "K/km", "K/m");

  m_precip_decay_rate = options::Real("-precip_decay_rate",
                                      "exp. decay rate for the surface mass balance,"
                                      "in 1/km",
                                      m_precip_decay_rate);
  m_precip_decay_rate = units::convert(m_sys, m_precip_decay_rate, "1 / km", "1 / m");

  m_precip_thresh_height = options::Real("-precip_thresh_height",
                                         "Threshold height for precipitation,"
                                         "km",
                                         m_precip_thresh_height);
  m_precip_thresh_height = units::convert(m_sys, m_precip_thresh_height, "km", "m");

  {
    ForcingOptions opt(*m_grid->ctx(), "atmosphere.index");

    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    m_log->message(2,
                  "  initializing %s, %s, %s and %s from forcing file %s...\n",
                  m_T0->get_name().c_str(),
                  m_T1->get_name().c_str(),
                  m_P0->get_name().c_str(),
                  m_P1->get_name().c_str(),
                  opt.filename.c_str());

    m_T0->init(opt.filename, 1, 0.0);
    m_T1->init(opt.filename, 1, 0.0);
    m_P0->init(opt.filename, 1, 0.0);
    m_P1->init(opt.filename, 1, 0.0);

    m_log->message(2,
                  "  reading %s data from forcing file %s...\n",
                  m_index->name().c_str(),
                  opt.filename.c_str());

    PIO nc(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);
    {
      m_index->read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
    }

    m_log->message(2,
                  "  reading surface elevation at t0 & t1 data from forcing file %s...\n",
                  opt.filename.c_str());

    m_h0.read(opt.filename, 0); // fails if not found!
    m_h1.read(opt.filename, 0);
  }
}



void IndexForcing::begin_pointwise_access_impl() const {
  m_T0->begin_access();
  m_T1->begin_access();
  m_P0->begin_access();
  m_P1->begin_access();
  m_h0.begin_access();
  m_h1.begin_access();
  m_surface->begin_access();
  m_precipitation.begin_access();
  m_air_temp.begin_access();
}

void IndexForcing::end_pointwise_access_impl() const {
  m_T0->end_access();
  m_T1->end_access();
  m_P0->end_access();
  m_P1->end_access();
  m_h0.end_access();
  m_h1.end_access();
  m_surface->end_access();
  m_precipitation.end_access();
  m_air_temp.end_access();
}

void IndexForcing::mean_annual_temp_impl(IceModelVec2S& result) const {
  result.copy_from(m_air_temp);
}

void IndexForcing::mean_precipitation_impl(IceModelVec2S& result) const {
  result.copy_from(m_precipitation);
}

void IndexForcing::update_impl(const Geometry &geometry ,double t, double dt) {
  const double t_index = m_grid->ctx()->time()->mod(t - m_index_reference_time, m_index_period),
               index = (*m_index)(t_index);

  m_T0->average(t, dt);
  m_T1->average(t, dt);
  m_P0->average(t, dt);
  m_P1->average(t, dt);

  m_surface = &(geometry.ice_surface_elevation);

  IceModelVec::AccessList list;
  list.add(m_air_temp);
  list.add(m_precipitation);
  list.add(*m_surface);
  list.add(*m_T0);
  list.add(*m_T1);
  list.add(*m_P0);
  list.add(*m_P1);
  list.add(m_h0);
  list.add(m_h1);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_air_temp(i, j) = compute_T_ij((*m_T0)(i, j), (*m_T1)(i, j), m_h0(i, j), m_h1(i, j), (*m_surface)(i, j), index);
    m_precipitation(i, j) = compute_P_ij((*m_P0)(i, j), (*m_P1)(i, j), m_h0(i, j), m_h1(i, j), (*m_surface)(i, j), index);
  }
}

void IndexForcing::define_model_state_impl(const PIO &output) const {
  m_precipitation.define(output);
  m_air_temp.define(output);
}

void IndexForcing::write_model_state_impl(const PIO &output) const {
  m_precipitation.write(output);
  m_air_temp.write(output);
}

void IndexForcing::init_timeseries_impl(const std::vector<double> &ts) const {
  int N = ts.size();
  m_ts_index.resize(N);
  m_ts_times.resize(N);

  for(unsigned int k = 0; k < N; k++){
    const double t = ts[k];
    m_ts_times[k] = t;
    m_ts_index[k] = (*m_index)(m_grid->ctx()->time()->mod(t - m_index_reference_time, m_index_period));
  }

  m_T0->init_interpolation(m_ts_times);
  m_T1->init_interpolation(m_ts_times);
  m_P0->init_interpolation(m_ts_times);
  m_P1->init_interpolation(m_ts_times);
}

void IndexForcing::temp_time_series_impl(int i, int j, std::vector<double>& result) const {
  std::vector<double> T0(m_ts_times.size()),
                      T1(m_ts_times.size());

  m_T0->interp(i, j, T0);
  m_T1->interp(i, j, T1);

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = compute_T_ij(T0[k], T1[k], m_h0(i, j), m_h1(i, j), (*m_surface)(i, j), m_ts_index[k]);
  }
}

void IndexForcing::precip_time_series_impl(int i, int j, std::vector<double>& result) const {
  std::vector<double> P0(m_ts_times.size()),
                      P1(m_ts_times.size());

  m_P0->interp(i, j, P0);
  m_P1->interp(i, j, P1);

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = compute_P_ij(P0[k], P1[k], m_h0(i, j), m_h1(i, j), (*m_surface)(i, j), m_ts_index[k]);
  }
}


double IndexForcing::compute_T_ij(double T0, double T1, double h0, double h1, double h, double index) const {
  const double T0_sl = applyLapseRateT(T0, h0, 0.0),
               T1_sl = applyLapseRateT(T1, h1, 0.0),
               T_sl = (T1_sl - T0_sl) * index + T0_sl,
               T = applyLapseRateT(T_sl, 0.0, h);

  return(T);
}

double IndexForcing::compute_P_ij(double P0, double P1, double h0, double h1, double h, double index) const {
  const double result = applyLapseRateP((P0 + (P1 - P0) * index), 0.0, h);
  return(std::max(0.0, result));
}

double IndexForcing::applyLapseRateT(double T, double h_ref, double h) const {
  const double result = T - m_temp_lapse_rate * (h - h_ref);
  return(result);
}

double IndexForcing::applyLapseRateP(double P, double h_ref, double h) const {
  (void) h_ref;
  const double result = P * exp(-1.0 * m_precip_decay_rate * std::max(0.0, (h - m_precip_thresh_height)));
  return(result);
}

} // end of namespace atmosphere
} // end of namespace pism
