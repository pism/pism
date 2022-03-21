// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2021, 2022 PISM Authors
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

#ifndef PISM_ATMOSPHERE_FRAC_P
#define PISM_ATMOSPHERE_FRAC_P

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {

class ScalarForcing;

namespace atmosphere {

class Frac_P : public AtmosphereModel {
public:
  Frac_P(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~Frac_P() = default;

private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  const array::Scalar& precipitation_impl() const;

  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;

  mutable std::vector<double> m_scaling_values;

  std::unique_ptr<ScalarForcing> m_1d_scaling;

  std::shared_ptr<array::Forcing> m_2d_scaling;

  array::Scalar::Ptr m_precipitation;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* PISM_ATMOSPHERE_FRAC_P */
