// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <cmath>

#include "base/part_grid_threshold_thickness.hh"
#include "base/util/Mask.hh"

namespace pism {


//! \file iMpartgrid.cc Methods implementing PIK option -part_grid [\ref Albrechtetal2011].

//! @brief Compute threshold thickness used when deciding if a
//! partially-filled cell should be considered 'full'.
double part_grid_threshold_thickness(StarStencil<int> M,
                                     StarStencil<double> H,
                                     StarStencil<double> h,
                                     double bed_elevation,
                                     double dx,
                                     bool reduce_frontal_thickness) {
  // get mean ice thickness and surface elevation over adjacent
  // icy cells
  double
    H_average   = 0.0,
    h_average   = 0.0,
    H_threshold = 0.0;
  int N = 0;

  const Direction dirs[] = {North, East, South, West};
  for (int n = 0; n < 4; ++n) {
    Direction direction = dirs[n];
    if (mask::icy(M[direction])) {
      H_average += H[direction];
      h_average += h[direction];
      N++;
    }
  }

  if (N == 0) {
    // If there are no "icy" neighbors, return the threshold thickness
    // of zero, forcing Href to be converted to H immediately.
    return 0.0;
  }

  H_average = H_average / N;
  h_average = h_average / N;

  if (bed_elevation + H_average > h_average) {
    H_threshold = h_average - bed_elevation;
  } else {
    H_threshold = H_average;
    // reduces the guess at the front
    if (reduce_frontal_thickness) {
      // FIXME: Magic numbers without references to the literature are bad.
      // for declining front C / Q0 according to analytical flowline profile in
      //   vandeveen with v0 = 300m year-1 and H0 = 600m
      const double
        H0         = 600.0,     // 600 m
        V0         = 300.0 / 3.15569259747e7, // 300 m year-1 (hard-wired for efficiency)
        mslope     = 2.4511e-18 * dx / (H0 * V0);
      H_threshold -= 0.8*mslope*pow(H_average, 5);
    }
  }

  // make sure that the returned threshold thickness is non-negative:
  return std::max(H_threshold, 0.0);
}

} // end of namespace pism
