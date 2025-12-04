/* Copyright (C) 2025 PISM Authors
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

#ifndef PISM_LONLATCALCULATOR_H
#define PISM_LONLATCALCULATOR_H

#include <array>

#include "pism/util/Proj.hh"

namespace pism {

/*!
 * Utility class converting `x,y` coordinates in a projection to a `lon,lat` pair.
 *
 * Requires the `PROJ` library.
 */
class LonLatCalculator {
public:
  LonLatCalculator(const std::string &proj_string)
      : m_coordinate_mapping(proj_string, "EPSG:4326") {
  }

  std::array<double, 2> lonlat(double x, double y) {
    PJ_COORD out = proj_trans(m_coordinate_mapping, PJ_FWD, proj_coord(x, y, 0, 0));

    // longitude: lambda
    // latitude: phi
    return { out.lp.lam, out.lp.phi };
  }

private:
  Proj m_coordinate_mapping;
};

} // namespace pism

#endif /* PISM_LONLATCALCULATOR_H */
