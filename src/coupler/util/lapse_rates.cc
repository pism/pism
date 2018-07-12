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

#include <cmath>                // fabs

#include "pism/util/iceModelVec.hh"

namespace pism {

void lapse_rate_scale(const IceModelVec2S &surface,
                      const IceModelVec2S &reference_surface,
                      double lapse_rate,
                      IceModelVec2S &result) {
  IceGrid::ConstPtr grid = result.grid();

  if (fabs(lapse_rate) < 1e-12) {
    return;
  }

  IceModelVec::AccessList list{&surface, &reference_surface, &result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    //as in paleo_precip
    const double correction = exp( -lapse_rate * (surface(i,j) - reference_surface(i,j)) );
    
    //as in paleo_precip, but taylor expansion
    //const double correction = 1.0 - lapse_rate * (surface(i,j) - reference_surface(i,j));

    //as in Pollard et al., 2012
    //not DeConto et al., 2016 ??
    //const double correction = pow(2.0, -lapse_rate/0.05/10.0 * (surface(i,j) - reference_surface(i,j)));

    result(i, j) *= correction;
  }
}

void lapse_rate_correction(const IceModelVec2S &surface,
                           const IceModelVec2S &reference_surface,
                           double lapse_rate,
                           IceModelVec2S &result) {
  IceGrid::ConstPtr grid = result.grid();

  if (fabs(lapse_rate) < 1e-12) {
    return;
  }

  IceModelVec::AccessList list{&surface, &reference_surface, &result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) -= lapse_rate * (surface(i,j) - reference_surface(i, j));
  }
}

} // end of namespace pism
