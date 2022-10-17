// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2021 PISM Authors
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

#ifndef PRECIPITATIONSCALING_H
#define PRECIPITATIONSCALING_H

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/ScalarForcing.hh"

namespace pism {

class ScalarForcing;

namespace atmosphere {

/** Scaling precipitation using air temperature offsets
 *
 * This class implements the adjustment used to capture the influence of changing air
 * temperature on precipitation. It uses scalar temperature offsets read from a file as an
 * input.
 *
 * @f[
 * P_{G}(x,y,t) = P_{G}(x,y,0)exp\left[\frac{0.169}{d}\left({\Delta}T(t)+ {\Delta}T_{SC}(t)\right)\right]
 * @f]

 * where @f$P_G(x,y,0)@f$ is the precipitation for the present
 * Greenland ice sheet, obtained from the atmosphere model used as an
 * input of this class. The time dependent precipitation is
 * @f$P_G(x,y,t)@f$, @f$d@f$ is the @f${\delta}_{18}@f$ conversion factor
 * of @f$2.4^{\circ}C/\frac{0}{00}@f$, and @f${\Delta}T_{SC}@f$ is a
 * correction term for the change of altitude of the central dome
 * during the Greenland ice sheet's evolution. The coefficient
 * " @f$ 0.169/d @f$ " corresponds to a 7.3% change of precipitation rate
 * for every @f$1^{\circ}C@f$ of temperature change (Huybrechts et al.
 * 2002).
 *
 * One possible scheme for @f${\Delta}T_{SC}(t)@f$ (used in PISM) is
 * to take it to be zero, which regards the height correction as
 * belonging to the set of uncertainties related to the conversion
 * between isotopic and temperature signals.
 */
class PrecipitationScaling : public AtmosphereModel {
public:
  PrecipitationScaling(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~PrecipitationScaling() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  void init_timeseries_impl(const std::vector<double> &ts) const;

  const array::Scalar& precipitation_impl() const;

  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;

protected:
  double m_exp_factor;
  std::unique_ptr<ScalarForcing> m_forcing;
  mutable std::vector<double> m_scaling_values;

  array::Scalar::Ptr m_precipitation;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* PRECIPITATIONSCALING_H */
