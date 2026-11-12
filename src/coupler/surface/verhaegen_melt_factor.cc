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

#include <algorithm>            // std::max
#include <cmath>                // std::pow

namespace pism {

double verhaegen_melt_factor(double debris_thickness) {

  // effective debris thickness (meters): the factor peaks here
  const double h_e = 0.015;
  // critical debris thickness (meters): the factor crosses 1 here
  const double h_c = 0.041;
  // lower bound corresponding to "no significant melt"
  const double f_min = 1e-3;

  auto &h = debris_thickness;

  // Clean ice. Written this way so that NaNs (and negative values produced by
  // interpolation) are treated as "no debris" instead of poisoning the melt rate.
  if (not (h > 0.0)) {
    return 1.0;
  }

  if (h <= h_e) {
    return 26.667 * h + 1.0;
  }

  if (h <= h_c) {
    return -16.0 * h + 1.64;
  }

  return std::max(f_min, 0.1061 * std::pow(h, -0.7205) - 0.07922);
}

}
