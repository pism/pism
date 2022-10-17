/* Copyright (C) 2014, 2015, 2016, 2017, 2018, 2021 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef PISM_WEATHER_STATION_HH
#define PISM_WEATHER_STATION_HH

#include "pism/coupler/AtmosphereModel.hh"

#include <memory>               // std::shared_ptr


namespace pism {

class ScalarForcing;

namespace atmosphere {

/** This class implements an atmosphere model corresponding to *one* weather station.
 *
 * It reads scalar (but time-dependent) near-surface air temperature
 * and precipitation data from a provided file. Resulting climate
 * fields are constant in time.
 *
 * This model should be used with a modifier such as `lapse_rate` to
 * create spatial variability.
 */
class WeatherStation : public AtmosphereModel {
public:
  WeatherStation(IceGrid::ConstPtr g);
  virtual ~WeatherStation() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const array::Scalar& precipitation_impl() const;
  const array::Scalar& air_temperature_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;
  void init_timeseries_impl(const std::vector<double> &ts) const;
  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &values) const;

  MaxTimestep max_timestep_impl(double t) const;
protected:

  std::shared_ptr<ScalarForcing> m_precipitation_timeseries;
  std::shared_ptr<ScalarForcing> m_air_temp_timeseries;

  mutable std::vector<double> m_precip_values, m_air_temp_values;

  array::Scalar::Ptr m_temperature;
  array::Scalar::Ptr m_precipitation;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* PISM_WEATHER_STATION_HH */
