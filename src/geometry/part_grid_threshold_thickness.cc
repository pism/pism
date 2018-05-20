// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/util/Mask.hh"

namespace pism {


//! \file part_grid_threshold_thickness.cc Methods implementing PIK option -part_grid [\ref Albrechtetal2011].

//! @brief Compute threshold thickness used when deciding if a
//! partially-filled cell should be considered 'full'.
double part_grid_threshold_thickness(StarStencil<int> M,
                                     StarStencil<double> H,
                                     StarStencil<double> h,
                                     double bed_elevation) {
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
  }

  // make sure that the returned threshold thickness is non-negative:
  return std::max(H_threshold, 0.0);
}

} // end of namespace pism
