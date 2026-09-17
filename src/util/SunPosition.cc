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

#include <cmath>

#include "pism/util/SunPosition.hh"

/*!
 * References:
 *
 * Sproul2007: A. B. Sproul, “Derivation of the solar geometric relationships using vector
 *             analysis,” Renewable Energy, vol. 32, no. 7, pp. 1187–1205, Jun. 2007, doi:
 *             10.1016/j.renene.2006.05.001.
 *
 * Zhang2021a: T. Zhang, P. W. Stackhouse, B. Macpherson, and J. C. Mikovitz, “A solar
 *             azimuth formula that renders circumstantial treatment unnecessary without
 *             compromising mathematical rigor: Mathematical setup, application and
 *             extension of a formula based on the subsolar point and atan2 function,”
 *             Renewable Energy, vol. 172, pp. 1333–1340, Jul. 2021, doi:
 *             10.1016/j.renene.2021.03.047.
 */

namespace pism {

static inline double clip(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

SunPosition::SunPosition(double declination) {
  m_sin_decl = std::sin(declination);
  m_cos_decl = std::cos(declination);
}

void SunPosition::set_hour_angles(const std::vector<double> &hour_angle) {
  auto N = hour_angle.size();
  m_cos_hour_angle.resize(N);
  m_sin_hour_angle.resize(N);
  for (int k = 0; k < N; ++k) {
    m_cos_hour_angle[k] = std::cos(hour_angle[k]);
    m_sin_hour_angle[k] = std::sin(hour_angle[k]);
  }
}

void SunPosition::set_latitude(double latitude_radians) {
  m_sin_lat = std::sin(latitude_radians);
  m_cos_lat = std::cos(latitude_radians);
}

void SunPosition::compute_at_set_hour_angle(int k, double &altitude, double &azimuth,
                                            double *solar_vector) const {
  compute_impl(m_cos_hour_angle[k], m_sin_hour_angle[k], altitude, azimuth, solar_vector);
}

void SunPosition::compute(double hour_angle, double &altitude, double &azimuth,
                          double *solar_vector) const {
  compute_impl(std::cos(hour_angle), std::sin(hour_angle), altitude, azimuth, solar_vector);
}

void SunPosition::compute_impl(double cos_hour_angle, double sin_hour_angle, double &altitude,
                               double &azimuth, double *solar_vector) const {

  // Vector pointing from the observer toward the center of the sun, equation (12) in
  // Sproul2007 or equations (6)--(8) in Zhang2021a:
  double S_e = - m_cos_decl * sin_hour_angle;
  double S_n = m_cos_lat * m_sin_decl - m_sin_lat * m_cos_decl * cos_hour_angle;
  double S_u = m_sin_lat * m_sin_decl + m_cos_lat * m_cos_decl * cos_hour_angle;

  if (solar_vector != nullptr) {
    solar_vector[0] = S_e;
    solar_vector[1] = S_n;
    solar_vector[2] = S_u;
  }

  // Altitude computed using equation (18) in Sproul2007 or equation (9) in Zhang2021a
  altitude = std::asin(clip(S_u, -1.0, 1.0));

  if (altitude >= 0.0) {
    // Equation (10) in Zhang2021a, corrected for the convention "north -> azimuth=0,
    // increase towards east of north" using text after equation (10).
    //
    // Note that x, y, z axes in Zhang2021a in fact point east, north, and up; see paragraph
    // about the derivation of equations (6)--(8).
    double A = std::atan2(S_e, S_n);

    // atan2(y,x) returns an angle in [-pi, pi]; here we adjust it to get azimuth in [0,
    // 2*pi):
    if (A < 0.0) {
      A += 2.0 * M_PI;
      // A tiny negative angle (e.g. at hour_angle = pi, where sin(hour_angle) is not
      // exactly zero in floating point) rounds to exactly 2*pi here; the documented range
      // is [0, 2*pi).
      if (A >= 2.0 * M_PI) {
        A = 0.0;
      }
    }
    azimuth = A;
  } else {
    azimuth = 0.0;
  }
}

} // end of namespace pism
