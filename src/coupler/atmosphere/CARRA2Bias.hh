// Copyright (C) 2026 PISM Authors
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

#ifndef PISM_ATMOSPHERE_CARRA2BIAS
#define PISM_ATMOSPHERE_CARRA2BIAS

#include <vector>

#include "pism/coupler/AtmosphereModel.hh"

namespace pism {
namespace atmosphere {

//! @brief Add an elevation- and season-dependent bias correction to the near-surface air
//! temperature.
//!
//! The bias is the difference between a free-air temperature carried down to the terrain
//! (a straight-line fit through CARRA2's pressure levels, evaluated at the surface) and
//! CARRA2's 2 m air temperature:
//!
//! @f[ b(z, t) = b_\infty + A \exp(-z / H)
//!     + (c_0 + c_1 z) \cos(2\pi t) + (s_0 + s_1 z) \sin(2\pi t), @f]
//!
//! where @f$ z @f$ is the surface elevation and @f$ t @f$ the fraction of the calendar
//! year. Below `z_min` and above `z_max` the elevation is held at the threshold, so the
//! bias is constant in elevation there (but keeps its seasonal cycle).
//!
//! This modifier returns @f$ T + b(z, t) @f$, where @f$ T @f$ is the air temperature of
//! the input model, a 2 m air temperature such as CARRA2's: @f$ b > 0 @f$ warms it toward
//! the free-air temperature at the surface.
class CARRA2Bias : public AtmosphereModel {
public:
  CARRA2Bias(std::shared_ptr<const Grid> g, std::shared_ptr<AtmosphereModel> in);
  virtual ~CARRA2Bias() = default;

  //! Bias (kelvin) at surface elevation `z` (meters), given cos(2 pi t) and sin(2 pi t).
  double bias(double z, double cos_t, double sin_t) const;

private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  const array::Scalar &air_temperature_impl() const;

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &values) const;

  //! Mean of cos(2 pi t) and sin(2 pi t) over the time step [t, t + dt].
  void mean_cycle(double t, double dt, double &cos_t, double &sin_t) const;

  // parameters, in kelvin, meters and kelvin per meter
  double m_background;
  double m_excess;
  double m_scale_height;
  double m_cos;
  double m_cos_gradient;
  double m_sin;
  double m_sin_gradient;
  double m_z_min;
  double m_z_max;

  // cos(2 pi t) and sin(2 pi t) at the times set by init_timeseries()
  mutable std::vector<double> m_cos_values;
  mutable std::vector<double> m_sin_values;

  std::shared_ptr<array::Scalar> m_temperature;

  // copy of the surface elevation, used by temp_time_series_impl()
  array::Scalar m_surface;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* PISM_ATMOSPHERE_CARRA2BIAS */
