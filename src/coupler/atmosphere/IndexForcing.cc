// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2023 PISM Authors
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

#include "IndexForcing.hh"

#include "pism/coupler/util/options.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace atmosphere {

IndexForcing::IndexForcing(IceGrid::ConstPtr grid)
  : AtmosphereModel(grid) {

  m_option = "atmosphere.index";

  m_ice_surface_elevation.create(m_grid, "m_ice_surface_elevation", WITHOUT_GHOSTS);
  m_ice_surface_elevation.set_attrs("internal", "ice surface elevation",
                                    "meters", "meters", "", 0);

  m_precip_exp_factor = m_config->get_number("atmosphere.precip_exponential_factor_for_temperature");

  m_temp_lapse_rate = m_config->get_number("atmosphere.elevation_change.temperature_lapse_rate",
                                           "K / m");

  m_precipitation = allocate_precipitation(grid);
  m_temperature = allocate_temperature(grid);

  m_h0.create(m_grid, "usurf_0", WITHOUT_GHOSTS);
  m_h0.set_attrs("climate_state", "surface elevation at t0",
                 "meters", "meters", "", 0);
  m_h0.set_time_independent(true);

  m_h1.create(m_grid, "usurf_1", WITHOUT_GHOSTS);
  m_h1.set_attrs("climate_state", "surface elevation at t1",
                 "meters", "meters", "", 0);
  m_h1.set_time_independent(true);

  m_index.reset(new ScalarForcing(grid->ctx(),
                                  m_option,
                                  "glac_index",
                                  "1",
                                  "1",
                                  "glacial index"));

  {
    auto filename = m_config->get_string(m_option + "_climate" + "_file");

    if (filename.empty()) {
      //If no extra file is specified, look in index file for climate data
      filename = m_config->get_string(m_option + ".file");
    }

    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_number("input.forcing.evaluations_per_year");
    const bool periodic = true;
    const InterpolationType interpolation_type = LINEAR_PERIODIC;

    File file(m_grid->com, filename, PISM_NETCDF3, PISM_READONLY);

    m_T0 = IceModelVec2T::ForcingField(m_grid,
                                       file,
                                       "airtemp_0",
                                       "", // no standard name
                                       buffer_size,
                                       evaluations_per_year,
                                       periodic,
                                       interpolation_type);
    m_T0->set_attrs("climate_forcing", "air temperature at t0",
                    "Kelvin", "Kelvin", "", 0);


    m_T1 = IceModelVec2T::ForcingField(m_grid,
                                       file,
                                       "airtemp_1",
                                       "", // no standard name
                                       buffer_size,
                                       evaluations_per_year,
                                       periodic,
                                       interpolation_type);
    m_T1->set_attrs("climate_forcing", "air temperature at t1",
                    "Kelvin", "Kelvin", "", 0);

    m_P0 = IceModelVec2T::ForcingField(m_grid,
                                       file,
                                       "precip_0",
                                       "", // no standard name
                                       buffer_size,
                                       evaluations_per_year,
                                       periodic,
                                       interpolation_type);
    m_P0->set_attrs("climate_forcing", "precipitation at t0",
                    "kg m-2 second-1", "kg m-2 second-1", "", 0);

    m_P1 = IceModelVec2T::ForcingField(m_grid,
                                       file,
                                       "precip_1",
                                       "", // no standard name
                                       buffer_size,
                                       evaluations_per_year,
                                       periodic,
                                       interpolation_type);
    m_P1->set_attrs("climate_forcing", "precipitation at t1",
                    "kg m-2 second-1", "kg m-2 second-1", "", 0);
  }
}

void IndexForcing::init_impl(const Geometry &geometry) {
  (void) geometry;

  auto filename = m_config->get_string(m_option + "_climate" + "_file");

  if (filename.empty()) {
    //If no extra file is specified, look in index file for climate data
    filename = m_config->get_string(m_option + ".file");
  }

  m_T0->init(filename, 1, 0.0);
  m_T1->init(filename, 1, 0.0);
  m_P0->init(filename, 1, 0.0);
  m_P1->init(filename, 1, 0.0);

  m_index->init();

  m_log->message(2,
                 "  reading surface elevation at t0 & t1 data from forcing file %s...\n",
                 filename.c_str());

  m_h0.regrid(filename, OPTIONAL, 0);
  m_h1.regrid(filename, OPTIONAL, 0);
}



void IndexForcing::begin_pointwise_access_impl() const {
  m_T0->begin_access();
  m_T1->begin_access();
  m_P0->begin_access();
  m_P1->begin_access();
  m_h0.begin_access();
  m_h1.begin_access();
  m_ice_surface_elevation.begin_access();
  m_precipitation->begin_access();
  m_temperature->begin_access();
}

void IndexForcing::end_pointwise_access_impl() const {
  m_T0->end_access();
  m_T1->end_access();
  m_P0->end_access();
  m_P1->end_access();
  m_h0.end_access();
  m_h1.end_access();
  m_ice_surface_elevation.end_access();
  m_precipitation->end_access();
  m_temperature->end_access();
}

const IceModelVec2S& IndexForcing::mean_precipitation_impl() const {
  return *m_precipitation;
}

const IceModelVec2S& IndexForcing::mean_annual_temp_impl() const {
  return *m_temperature;
}

void IndexForcing::update_impl(const Geometry &geometry ,double t, double dt) {
  m_index->update(t, dt);

  m_T0->average(t, dt);
  m_T1->average(t, dt);
  m_P0->average(t, dt);
  m_P1->average(t, dt);
  m_ice_surface_elevation.copy_from(geometry.ice_surface_elevation);

  IceModelVec::AccessList list{ m_temperature.get(), m_precipitation.get(),
                                m_T0.get(), m_T1.get(), m_P0.get(), m_P1.get(),
                                &m_ice_surface_elevation, &m_h0, &m_h1 };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*m_temperature)(i, j)   = compute_T_ij((*m_T0)(i, j), (*m_T1)(i, j),
                                            m_h0(i, j), m_h1(i, j),
                                            m_ice_surface_elevation(i, j),
                                            m_index->value());
    (*m_precipitation)(i, j) = compute_P_ij((*m_P0)(i, j), (*m_P1)(i, j),
                                            m_h0(i, j), m_h1(i, j),
                                            m_ice_surface_elevation(i, j),
                                            m_index->value());
  }
}

void IndexForcing::define_model_state_impl(const File &output) const {
  m_precipitation->define(output);
  m_temperature->define(output);
}

void IndexForcing::write_model_state_impl(const File &output) const {
  m_precipitation->write(output);
  m_temperature->write(output);
}

void IndexForcing::init_timeseries_impl(const std::vector<double> &ts) const {
  int N = ts.size();
  m_ts_index.resize(N);
  m_ts_times.resize(N);

  for(int k = 0; k < N; k++){
    const double t = ts[k];
    m_ts_times[k] = t;
    m_ts_index[k] = m_index->value(t);
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

  const double h0_ij = m_h0(i, j),
               h1_ij = m_h1(i, j),
               ice_surface_elevation_ij = m_ice_surface_elevation(i, j);

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = compute_T_ij(T0[k], T1[k], h0_ij, h1_ij, ice_surface_elevation_ij, m_ts_index[k]);
  }
}

void IndexForcing::precip_time_series_impl(int i, int j, std::vector<double>& result) const {
  std::vector<double> P0(m_ts_times.size()),
                      P1(m_ts_times.size());

  m_P0->interp(i, j, P0);
  m_P1->interp(i, j, P1);

  const double h0_ij = m_h0(i, j),
               h1_ij = m_h1(i, j),
               ice_surface_elevation_ij = m_ice_surface_elevation(i, j);

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = compute_P_ij(P0[k], P1[k], h0_ij, h1_ij, ice_surface_elevation_ij, m_ts_index[k]);
  }
}


double IndexForcing::compute_T_ij(double T0, double T1, double h0, double h1, double h, double index) const {
  const double T0_sl = applyLapseRateT(T0, h0, 0.0),
               T1_sl = applyLapseRateT(T1, h1, 0.0),
               T_sl  = (T1_sl - T0_sl) * index + T0_sl,
               T     = applyLapseRateT(T_sl, 0.0, h);

  return T;
}

double IndexForcing::compute_P_ij(double P0, double P1, double h0, double h1, double h, double index) const {
  const double P0_sl = applyLapseRateP(P0, h0, 0.0),
               P1_sl = applyLapseRateP(P1, h1, 0.0),
               P_sl  = (P1_sl - P0_sl) * index + P0_sl,
               P     = applyLapseRateP(P_sl, 0.0, h);

  return std::max(0.0, P);
}

double IndexForcing::applyLapseRateT(double T, double h_ref, double h) const {
  const double result = T - m_temp_lapse_rate * (h - h_ref);
  return(result);
}

double IndexForcing::applyLapseRateP(double P, double h_ref, double h) const {
  const double result = P * exp(-1.0  * m_precip_exp_factor * m_temp_lapse_rate * (h - h_ref));
  return result;
}

} // end of namespace atmosphere
} // end of namespace pism
