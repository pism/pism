/* Copyright (C) 2026 PISM Authors
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

#include <vector>

namespace pism {

/*!
 * Compute sun position (altitude and azimuth in radians) given declination, latitude, and
 * hour angle.
 *
 * Written in a way allowing us to avoid unnecessarily re-computing sin(declination),
 * cos(declination) more than once per time step (these depend only on time) and
 * sin(latitude) and cos(latitude) (these do not depend on hour angle or time, but we
 * re-compute them once per time step anyway).
 *
 * Call sequence:
 *
 * 1. Create an instance:
 *
 * SunPosition sun_position(declination);
 *
 * 2. Set latitude:
 *
 * sun_position.set_latitude(latitude);
 *
 * 3. Set hour angles to pre-compute sin() and cos():
 *
 * sun_position.set_hour_angles(angles);
 *
 * 3. Compute sun position in the sky given an hour angle index:
 *
 * double altitude, azimuth;
 * double solar_vector[3];
 * sun_position.compute_at_set_hour_angle(k, altitude, azimuth, solar_vector);
 *
 * Azimuth is set to zero if altitude < 0.
 */
class SunPosition {
public:
  SunPosition(double declination);

  void set_latitude(double latitude_radians);

  void set_hour_angles(const std::vector<double> &hour_angle);
  void compute_at_set_hour_angle(int k, double &altitude, double &azimuth,
                                 double *solar_vector = nullptr) const;

  void compute(double hour_angle, double &altitude, double &azimuth,
               double *solar_vector = nullptr) const;

private:
  double m_sin_decl, m_cos_decl;
  double m_sin_lat, m_cos_lat;

  void compute_impl(double cos_hour_angle, double sin_hour_angle, double &altitude,
                    double &azimuth, double *solar_vector) const;

  std::vector<double> m_sin_hour_angle;
  std::vector<double> m_cos_hour_angle;
};

} // end of namespace pism
