// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {
namespace atmosphere {

//! \brief Reads and uses air_temp and precipitation anomalies from a file.
class Anomaly : public AtmosphereModel
{
public:
  Anomaly(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~Anomaly();

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& mean_precipitation_impl() const;
  const IceModelVec2S& mean_annual_temp_impl() const;

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;
  void temp_time_series_impl(int i, int j, std::vector<double> &values) const;
  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;
protected:
  mutable std::vector<double> m_mass_flux_anomaly, m_temp_anomaly;

  IceModelVec2T::Ptr m_air_temp_anomaly;
  IceModelVec2T::Ptr m_precipitation_anomaly;

  IceModelVec2S::Ptr m_precipitation;
  IceModelVec2S::Ptr m_temperature;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAANOMALY_H_ */
