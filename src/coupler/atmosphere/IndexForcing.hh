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

#include "pism/coupler/AtmosphereModel.hh"

#include "pism/coupler/util/ScalarForcing.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {
namespace atmosphere {

class IndexForcing : public AtmosphereModel {

public:
  IndexForcing(IceGrid::ConstPtr g);

protected:
  std::unique_ptr<ScalarForcing> m_index;
  unsigned int m_index_period;
  double m_index_reference_time;
  double m_precip_exp_factor;
  double m_temp_lapse_rate;
  mutable std::vector<double> m_ts_index;
  IceModelVec2S m_h0,
                m_h1;
  IceModelVec2S m_ice_surface_elevation;
  IceModelVec2T::Ptr m_T0,
                     m_T1,
                     m_P0,
                     m_P1;

  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);
  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  const IceModelVec2S& mean_precipitation_impl() const;
  const IceModelVec2S& mean_annual_temp_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;
  void init_timeseries_impl(const std::vector<double> &ts) const;
  void precip_time_series_impl(int i, int j, std::vector<double> &result) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &result) const;

private:
  IceModelVec2S::Ptr m_precipitation, m_temperature;
  double applyLapseRateT(double T, double h_ref, double h) const;
  double applyLapseRateP(double P, double h_ref, double h) const;
  double compute_T_ij(double T0, double T1, double h0, double h1, double h, double index) const;
  double compute_P_ij(double P0, double P1, double h0, double h1, double h, double index) const;
  std::string m_option;

};

} // end of namespace atmosphere
} // end of namespace pism
