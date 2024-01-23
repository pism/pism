// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023 PISM Authors
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

#ifndef _CLIMATEINDEX_H_
#define _CLIMATEINDEX_H_

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/coupler/util/ClimateIndexWeights.hh"

namespace pism {

class ClimateIndexWeights;
class ScalarForcing;

namespace atmosphere {

class ClimateIndex : public AtmosphereModel {
public:
  ClimateIndex(std::shared_ptr<const Grid> g);
  virtual ~ClimateIndex() = default;

  virtual const array::Scalar& mean_summer_temp() const;

protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual const array::Scalar& precipitation_impl() const;
  virtual const array::Scalar& air_temperature_impl() const;

  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;

  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void temp_time_series_impl(int i, int j, std::vector<double> &result) const;
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &result) const;

  virtual DiagnosticList diagnostics_impl() const;
protected:
  std::unique_ptr<ScalarForcing> m_A; // amplitude scaling
  std::unique_ptr<ClimateIndexWeights> m_climate_index;

  bool use_cos, use_1X;
  double m_w0, m_w1, m_w1X;
  double m_snow_temp_summer_day;
  double m_preciplinfactor;
  bool use_precip_scaling;
  std::string m_reference, precip_scaling_file;
  mutable std::vector<double> m_cosine_cycle;

  array::Scalar m_air_temp_annual, m_air_temp_annual_ref, m_air_temp_summer, m_air_temp_summer_ref;
  array::Scalar m_air_temp_anomaly_annual_0, m_air_temp_anomaly_annual_1,m_air_temp_anomaly_annual_1X;
  array::Scalar m_air_temp_anomaly_summer_0, m_air_temp_anomaly_summer_1, m_air_temp_anomaly_summer_1X;

  array::Scalar m_precipitation, m_precipitation_ref;
  array::Scalar m_precipitation_anomaly_0, m_precipitation_anomaly_1, m_precipitation_anomaly_1X;
  array::Scalar m_spatial_precip_scaling;

};

} // end of namespace atmosphere
} // end of namespace pism


#endif /* _CLIMATEINDEX_H_ */
