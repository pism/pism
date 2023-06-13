/* Copyright (C) 2018, 2022 PISM Authors
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

#include <cmath>                // fabs

#include "pism/util/array/Scalar.hh"

namespace pism {

void lapse_rate_correction(const array::Scalar &surface,
                           const array::Scalar &reference_surface,
                           double lapse_rate,
                           array::Scalar &result) {
  IceGrid::ConstPtr grid = result.grid();

  if (fabs(lapse_rate) < 1e-12) {
    return;
  }

  array::AccessScope list{&surface, &reference_surface, &result};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) -= lapse_rate * (surface(i,j) - reference_surface(i, j));
  }
}

} // end of namespace pism
