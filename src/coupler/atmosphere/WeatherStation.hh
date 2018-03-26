/* Copyright (C) 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PAWEATHERSTATION_H_
#define _PAWEATHERSTATION_H_

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/Timeseries.hh"

namespace pism {
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
  virtual ~WeatherStation();

protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual const IceModelVec2S& mean_precipitation_impl() const;
  virtual const IceModelVec2S& mean_annual_temp_impl() const;

  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;
  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &values) const;
  virtual void temp_time_series_impl(int i, int j, std::vector<double> &values) const;

  virtual MaxTimestep max_timestep_impl(double t) const;
protected:
  Timeseries m_precipitation_timeseries, m_air_temp_timeseries;
  mutable std::vector<double> m_precip_values, m_air_temp_values;

  IceModelVec2S::Ptr m_temperature;
  IceModelVec2S::Ptr m_precipitation;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAWEATHERSTATION_H_ */
