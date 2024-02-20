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

#ifndef PISM_ATMOSPHERE_CLIMATEINDEX_H
#define PISM_ATMOSPHERE_CLIMATEINDEX_H

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/coupler/util/ClimateIndexWeights.hh"

namespace pism {

class ScalarForcing;

namespace atmosphere {

class ClimateIndex : public AtmosphereModel {
public:
  ClimateIndex(std::shared_ptr<const Grid> g);
  virtual ~ClimateIndex() = default;

  virtual const array::Scalar& mean_summer_temp() const;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const array::Scalar& precipitation_impl() const;
  const array::Scalar& air_temperature_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &result) const;
  void precip_time_series_impl(int i, int j, std::vector<double> &result) const;

  DiagnosticList diagnostics_impl() const;

  std::unique_ptr<ScalarForcing> m_A; // amplitude scaling
  std::unique_ptr<ClimateIndexWeights> m_climate_index;

  bool m_use_cos, m_use_precip_cos, m_use_1X;
  double m_midsummer_year_fraction;

  bool m_use_precip_scaling;
  bool m_spatially_variable_scaling;
  double m_preciplinfactor;

  mutable std::vector<double> m_cosine_cycle;

  array::Scalar m_air_temp_annual;
  array::Scalar m_air_temp_annual_ref;
  array::Scalar m_air_temp_summer;
  array::Scalar m_air_temp_summer_ref;
  array::Scalar m_air_temp_anomaly_annual_0;
  array::Scalar m_air_temp_anomaly_annual_1;
  array::Scalar m_air_temp_anomaly_annual_1X;
  array::Scalar m_air_temp_anomaly_summer_0;
  array::Scalar m_air_temp_anomaly_summer_1;
  array::Scalar m_air_temp_anomaly_summer_1X;

  array::Scalar m_precipitation;
  array::Scalar m_precipitation_ref;
  array::Scalar m_precipitation_anomaly_0;
  array::Scalar m_precipitation_anomaly_1;
  array::Scalar m_precipitation_anomaly_1X;
  array::Scalar m_spatial_precip_scaling;
};

} // end of namespace atmosphere
} // end of namespace pism


#endif /* PISM_ATMOSPHERE_CLIMATEINDEX_H */
