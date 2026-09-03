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

#ifndef PISM_MPDATA_HELPERS_H
#define PISM_MPDATA_HELPERS_H

#include <algorithm>            // std::min, std::max

namespace pism {
namespace mpdata {

//! positive part
inline double pp(double x) {
  return std::max(x, 0.0);
}

//! negative part
inline double np(double x) {
  return std::min(x, 0.0);
}

//! Flux across a face with velocity `u` using first-order upwinding (`x` on the side the
//! velocity points away from, `x_n` on the other side).
inline double upwind(double x, double x_n, double u) {
  return u * (u >= 0.0 ? x : x_n);
}

/*!
 * Anti-diffusive ("corrective") velocity of Smolarkiewicz (1983), equation 13, for a
 * face between two cells with values `x_c` and `x_n` separated by `d`, given the
 * upwind-pass velocity `u` across the face.
 */
inline double corrective_velocity(double dt, double u, double x_c, double x_n, double d,
                                  double eps) {
  return (std::abs(u) * d - dt * u * u) * (x_n - x_c) / ((x_n + x_c + eps) * d);
}

/*!
 * Cross term of the anti-diffusive velocity: `-0.5 dt u v_bar (dx/dy) / x_bar`, where the
 * derivative in the transverse direction is approximated by the difference of the sums
 * `s_plus - s_minus` of the values on the two sides, `s_total` is their sum, and `d` is
 * the transverse spacing.
 */
inline double cross_term(double dt, double u, double v_bar, double s_plus, double s_minus,
                         double s_total, double d, double eps) {
  return -0.5 * dt * u * v_bar * (s_plus - s_minus) / ((s_total + eps) * d);
}

} // end of namespace mpdata
} // end of namespace pism

#endif /* PISM_MPDATA_HELPERS_H */
