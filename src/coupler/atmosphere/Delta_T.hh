// Copyright (C) 2011--2022 PISM Authors
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

#ifndef PISM_ATMOSPHERE_DELTA_T
#define PISM_ATMOSPHERE_DELTA_T

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {

class ScalarForcing;

namespace atmosphere {

class Delta_T : public AtmosphereModel {
public:
  Delta_T(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~Delta_T() = default;

private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  const array::Scalar& air_temperature_impl() const;

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &values) const;

  mutable std::vector<double> m_offset_values;

  std::unique_ptr<ScalarForcing> m_1d_offsets;

  std::shared_ptr<IceModelVec2T> m_2d_offsets;

  array::Scalar::Ptr m_temperature;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* PISM_ATMOSPHERE_DELTA_T */
