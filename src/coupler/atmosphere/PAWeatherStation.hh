/* Copyright (C) 2014 PISM Authors
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

#include "PISMAtmosphere.hh"
#include "Timeseries.hh"

namespace pism {

/** This class implements an atmosphere model corresponding to *one* weather station.
 *
 * It reads scalar (but time-dependent) near-surface air temperature
 * and precipitation data from a provided file. Resulting climate
 * fields are constant in time.
 *
 * This model should be used with a modifier such as `lapse_rate` to
 * create spatial variability.
 */
class PAWeatherStation : public AtmosphereModel {
public:
  PAWeatherStation(IceGrid &g);
  virtual ~PAWeatherStation();

  virtual void init();

  virtual void update(double t, double dt);
  virtual void mean_precipitation(IceModelVec2S &result);
  virtual void mean_annual_temp(IceModelVec2S &result);

  virtual void begin_pointwise_access();
  virtual void end_pointwise_access();
  virtual void init_timeseries(const std::vector<double> &ts);
  virtual void precip_time_series(int i, int j, std::vector<double> &values);
  virtual void temp_time_series(int i, int j, std::vector<double> &values);
  virtual void temp_snapshot(IceModelVec2S &result);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);

  virtual void write_variables(const std::set<std::string> &vars, const PIO& nc);

protected:
  Timeseries m_precipitation, m_air_temperature;
  std::vector<double> m_precip_values, m_air_temp_values;

  NCSpatialVariable m_precip_metadata, m_air_temp_metadata;
};

}      // end of namespace pism

#endif /* _PAWEATHERSTATION_H_ */
