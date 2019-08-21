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

#ifndef _PAYEARLYCYCLE_H_
#define _PAYEARLYCYCLE_H_

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace atmosphere {

//! A class containing an incomplete implementation of an atmosphere model
//! based on a temperature parameterization using mean annual and mean
//! summer temperatures and a cosine yearly cycle. Uses a stored
//! (constant in time) precipitation field.
class YearlyCycle : public AtmosphereModel {
public:
  YearlyCycle(IceGrid::ConstPtr g);
  virtual ~YearlyCycle();

  virtual const IceModelVec2S& mean_summer_temp() const;

protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual void init_impl(const Geometry &geometry);

  virtual const IceModelVec2S& mean_precipitation_impl() const;
  virtual const IceModelVec2S& mean_annual_temp_impl() const;

  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;

  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual void temp_time_series_impl(int i, int j, std::vector<double> &result) const;
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &result) const;

  virtual void update_impl(const Geometry &geometry, double t, double dt) = 0;

  virtual DiagnosticList diagnostics_impl() const;
protected:
  void init_internal(const std::string &input_filename, bool regrid,
                     unsigned int start);

  double m_snow_temp_summer_day;
  std::string m_reference;
  IceModelVec2S m_air_temp_mean_annual, m_air_temp_mean_summer, m_precipitation;
  mutable std::vector<double> m_ts_times;
  mutable std::vector<double> m_cosine_cycle;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAYEARLYCYCLE_H_ */
