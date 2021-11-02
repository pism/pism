// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2021 PISM Authors
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

#ifndef _PAELEVATIONCHANGE_H_
#define _PAELEVATIONCHANGE_H_

#include "pism/coupler/AtmosphereModel.hh"

#include "pism/util/iceModelVec2T.hh"

namespace pism {
namespace atmosphere {

class ElevationChange : public AtmosphereModel
{
public:
  ElevationChange(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~ElevationChange() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& mean_precipitation_impl() const;
  const IceModelVec2S& mean_annual_temp_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void precip_time_series_impl(int i, int j, std::vector<double> &result) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &result) const;

protected:
  enum Method {SCALE, SHIFT};

  Method m_precip_method;

  // parameter used when m_precip_method == SHIFT
  double m_precip_lapse_rate;

  // parameters used when m_precip_method == SCALE
  double m_precip_temp_lapse_rate;
  double m_precip_exp_factor;

  // surface temperature lapse rate
  double m_temp_lapse_rate;

  std::shared_ptr<IceModelVec2T> m_reference_surface;

  IceModelVec2S::Ptr m_precipitation;
  IceModelVec2S::Ptr m_temperature;
  IceModelVec2S m_surface;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAELEVATIONCHANGE_H_ */
