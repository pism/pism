/* Copyright (C) 2018 PISM Authors
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

#include "pism/coupler/AtmosphereModel.hh"

namespace pism {
namespace atmosphere {

/*!
 * Test atmosphere model that returns air temperature and precipitation fields that are
 * constant in time and space (for testing only).
 */
class Uniform : public AtmosphereModel {
public:
  Uniform(IceGrid::ConstPtr g);
private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& mean_precipitation_impl() const;
  const IceModelVec2S& mean_annual_temp_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  void temp_time_series_impl(int i, int j, std::vector<double> &values) const;
  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;

private:
  IceModelVec2S::Ptr m_precipitation, m_temperature;
};

} // end of namespace atmosphere
} // end of namespace pism
