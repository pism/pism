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

void SunPosition::compute_at_set_hour_angle(int k, double &altitude, double &azimuth) const {
  compute_impl(m_cos_hour_angle[k], m_sin_hour_angle[k], altitude, azimuth);
}

void SunPosition::compute(double hour_angle, double &altitude, double &azimuth) const {
  compute_impl(std::cos(hour_angle), std::sin(hour_angle), altitude, azimuth);
}

// Standard topocentric solar geometry (textbook spherical astronomy). The ENU sun-vector
// convention (E = cos(alt) sin(az), N = cos(alt) cos(az), U = sin(alt)) matches solshade's
// solar.py; the altitude/azimuth formulas themselves are standard.
void SunPosition::compute_impl(double cos_hour_angle, double sin_hour_angle, double &altitude,
                               double &azimuth) const {
  double sin_alt =
      clip(m_sin_lat * m_sin_decl + m_cos_lat * m_cos_decl * cos_hour_angle, -1.0, 1.0);

  altitude = std::asin(sin_alt);

  // Azimuth is irrelevant if the sun is below the horizon.
  if (altitude < 0.0) {
    azimuth = 0.0;
    return;
  }

  double cos_altitude = std::cos(altitude);

  // Degenerate geometry: sun at the zenith, or observer at a geographic pole. Azimuth is
  // undefined; return 0 (irrelevant for the cosine projection at the zenith, and PISM
  // domains are not located exactly at a pole).
  if (cos_altitude < 1e-8 || m_cos_lat < 1e-8) {
    azimuth = 0.0;
    return;
  }

  double sinA = -m_cos_decl * sin_hour_angle / cos_altitude;
  double cosA = (m_sin_decl - m_sin_lat * sin_alt) / (m_cos_lat * cos_altitude);

  double A = std::atan2(sinA, cosA); // clockwise from north, in (-pi, pi]

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
}

} // end of namespace pism
