// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#ifndef _PAINDEXFORCING_H_
#define _PAINDEXFORCING_H_

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {
namespace atmosphere {

class IndexForcing : public AtmosphereModel {

public:
  IndexForcing(IceGrid::ConstPtr g);
  ~IndexForcing();

protected:
  unsigned int m_index_period;
  double m_index_reference_time, m_temp_lapse_rate, m_precip_decay_rate,
         m_precip_thresh_height;
  mutable std::vector<double> m_ts_index;
  Timeseries *m_index;
  IceModelVec2S m_precipitation, m_air_temp, m_h0, m_h1;
  const IceModelVec2S *m_surface;
  IceModelVec2T::Ptr m_T0, m_T1, m_P0, m_P1;

  virtual void init_impl(const Geometry &geometry);
  virtual void mean_precipitation_impl(IceModelVec2S &result) const;
  virtual void mean_annual_temp_impl(IceModelVec2S &result) const;
  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;
  virtual void temp_time_series_impl(int i, int j, std::vector<double> &values) const;
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &values) const;
  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;
  virtual void update_impl(const Geometry &geometry, double t, double dt);

private:
  double applyLapseRateT(double T, double h_ref, double h) const;
  double applyLapseRateP(double P, double h_ref, double h) const;
  double compute_T_ij(double T0, double T1, double h0, double h1, double h, double index) const;
  double compute_P_ij(double P0, double P1, double h0, double h1, double h, double index) const;

};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAINDEXFORCING_H_ */
