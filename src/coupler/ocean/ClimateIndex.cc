// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023, 2024 PISM Authors
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

#include "ClimateIndex.hh"
#include "pism/util/Grid.hh"

namespace pism {
namespace ocean {

ClimateIndex::ClimateIndex(std::shared_ptr<const Grid> grid)
    : CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),

      m_theta_ref(m_grid, "theta_ocean_ref"),
      m_salinity_ref(m_grid, "salinity_ocean_ref"),

      m_theta_anomaly_0(m_grid, "theta_ocean_anomaly_0"),
      m_salinity_anomaly_0(m_grid, "salinity_ocean_anomaly_0"),

      m_theta_anomaly_1(m_grid, "theta_ocean_anomaly_1"),
      m_salinity_anomaly_1(m_grid, "salinity_ocean_anomaly_1"),

      m_theta_anomaly_1X(m_grid, "theta_ocean_anomaly_1X"),
      m_salinity_anomaly_1X(m_grid, "salinity_ocean_anomaly_1X") {

  auto climate_index_file = m_config->get_string("climate_index.file");
  if (climate_index_file.empty()) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION, "'Climate Index Weight' file `climate_index.file` cannot be empty");
  }

  m_climate_index.reset(new ClimateIndexWeights(*grid->ctx()));

  m_theta_ref.metadata(0)
      .long_name("potential temperature of the adjacent ocean")
      .units("Kelvin")
      .set_time_independent(true);

  // Paleo time slice temperature data annual
  m_theta_anomaly_0.metadata(0)
      .long_name("absolute potential temperature anomaly of the adjacent ocean")
      .units("Kelvin")
      .set_time_independent(true);

  m_theta_anomaly_1.metadata(0)
      .long_name("potential temperature of the adjacent ocean")
      .units("Kelvin")
      .set_time_independent(true);

  m_theta_anomaly_1X.metadata(0)
      .long_name("potential temperature of the adjacent ocean")
      .units("Kelvin")
      .set_time_independent(true);

  m_salinity_ref.metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg")
      .set_time_independent(true);

  // Paleo time slice temperature data annual
  m_salinity_anomaly_0.metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg")
      .set_time_independent(true);

  m_salinity_anomaly_1.metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg")
      .set_time_independent(true);

  m_salinity_anomaly_1X.metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg")
      .set_time_independent(true);

  m_use_1X = m_config->get_flag("climate_index.super_interglacial.use");
}

void ClimateIndex::init_forcing() {
  m_log->message(2, "**** Initializing the 'Climate Index' ocean model...\n");

  auto input_file = m_config->get_string("ocean.climate_index.climate_snapshots.file");
  if (input_file.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "An ocean snapshots input file name\n"
                       "ocean.climate_index.climate_snapshots.file cannot be empty");
  }
  m_log->message(2,
                 "  Reading salinity and theta fields from '%s'...\n",
                 input_file.c_str());

  File input(m_grid->com, input_file, io::PISM_GUESS, io::PISM_READONLY);
  auto None = io::Default::Nil();

  // Reference fields
  m_theta_ref.regrid(input, None);
  m_salinity_ref.regrid(input, None);

  // Annual anomaly for Paleo time slices 0=Glacial, 1=Interglacial, 1X= Super InterGlacial e.g. mPWP
  m_theta_anomaly_0.regrid(input, None);

    m_log->message(2,
            " m_theta_ocean_0 loaded ");
  m_theta_anomaly_1.regrid(input, None);

  m_log->message(2,
            " m_theta_ocean_1 loaded ");

  m_salinity_anomaly_0.regrid(input, None);
    m_log->message(2,
            " m_salinity_ocean_1 loaded ");
  m_salinity_anomaly_1.regrid(input, None);

  if (m_use_1X) {
    m_theta_anomaly_1X.regrid(input, None);
    m_salinity_anomaly_1X.regrid(input, None);
  }
}

void ClimateIndex::update_forcing(double t, double dt, array::Scalar &theta,
                                  array::Scalar &salinity) {

  auto weights = m_climate_index->update_weights(t, dt);
  double w_0   = weights[0];
  double w_1   = weights[1];
  double w_1X  = m_use_1X ? weights[2] : 0.0;

  m_log->message(3, "**** Updated weights in ocean: w0 = '%f', w1 = '%f', w1X = '%f' ****\n",
                 w_0, w_1, w_1X);

  const auto &T_ref = m_theta_ref;
  const auto &dT_0  = m_theta_anomaly_0;
  const auto &dT_1  = m_theta_anomaly_1;
  const auto &dT_1X = m_theta_anomaly_1X;

  const auto &S_ref = m_salinity_ref;
  const auto &dS_0  = m_salinity_anomaly_0;
  const auto &dS_1  = m_salinity_anomaly_1;
  const auto &dS_1X = m_salinity_anomaly_1X;

  array::AccessScope scope{ &theta, &salinity, &T_ref, &dT_0,  &dT_1,
                            &S_ref, &dS_0,     &dS_1,  &dT_1X, &dS_1X };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    theta(i, j) =
        T_ref(i, j) + w_0 * dT_0(i, j) + w_1 * dT_1(i, j) + w_1X * (dT_1X(i, j) - dT_1(i, j));
    salinity(i, j) =
        S_ref(i, j) + w_0 * dS_0(i, j) + w_1 * dS_1(i, j) + w_1X * (dS_1X(i, j) - dS_1(i, j));
  }
}

} // end of namespace ocean
} // end of namespace pism
