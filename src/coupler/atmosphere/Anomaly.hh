// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#ifndef _PAANOMALY_H_
#define _PAANOMALY_H_

#include "pism/coupler/util/PGivenClimate.hh"
#include "Modifier.hh"

namespace pism {
namespace atmosphere {

//! \brief Reads and uses air_temp and precipitation anomalies from a file.
class Anomaly : public PGivenClimate<PAModifier,AtmosphereModel>
{
public:
  Anomaly(IceGrid::ConstPtr g, AtmosphereModel* in);
  virtual ~Anomaly();

protected:
  virtual void init_impl();
  virtual void update_impl(double my_t, double my_dt);

  virtual void mean_precipitation_impl(IceModelVec2S &result) const;
  virtual void mean_annual_temp_impl(IceModelVec2S &result) const;

  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;
  virtual void temp_time_series_impl(int i, int j, std::vector<double> &values) const;
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &values) const;
protected:
  bool m_modify_precip;
  std::vector<double> m_ts_mod, m_ts_values;
  IceModelVec2T *m_air_temp_anomaly, *m_precipitation_anomaly;
  mutable std::vector<double> m_mass_flux_anomaly, m_temp_anomaly;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAANOMALY_H_ */
