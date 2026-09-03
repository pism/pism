// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#include <algorithm>            // std::max
#include <cmath>                // std::sqrt

#include "pism/coupler/debris/terminus_kernels.hh"

namespace pism {
namespace debris {
namespace terminus {

std::vector<Offset> upstream_chain(const std::function<double(int, int)> &surface,
                                   const std::function<bool(int, int)> &icy,
                                   int i, int j, int n, double dx, double dy) {
  std::vector<Offset> result;

  int ci = i, cj = j;
  for (int step = 0; step < n; ++step) {
    const double h_c = surface(ci, cj);
    double slope_max = 0.0;
    int best_i = ci, best_j = cj;

    for (int dj = -1; dj <= 1; ++dj) {
      for (int di = -1; di <= 1; ++di) {
        if (di == 0 and dj == 0) {
          continue;
        }
        const int ni = ci + di, nj = cj + dj;
        if (not icy(ni, nj)) {
          continue;
        }
        const double distance = std::sqrt(di * di * dx * dx + dj * dj * dy * dy);
        const double slope = (surface(ni, nj) - h_c) / distance;
        if (slope > slope_max) {
          slope_max = slope;
          best_i    = ni;
          best_j    = nj;
        }
      }
    }

    if (best_i == ci and best_j == cj) {
      // local maximum: stop
      break;
    }

    ci = best_i;
    cj = best_j;
    result.push_back({ ci - i, cj - j });
  }

  return result;
}

double foreland_slope(const std::function<double(int, int)> &surface,
                      const std::function<bool(int, int)> &ice_free,
                      int i, int j, double dx, double dy) {
  const double h = surface(i, j);

  double s_x = 0.0, s_y = 0.0;

  for (int di : { -1, 1 }) {
    if (ice_free(i + di, j)) {
      s_x = std::max(s_x, (h - surface(i + di, j)) / dx);
    }
  }
  for (int dj : { -1, 1 }) {
    if (ice_free(i, j + dj)) {
      s_y = std::max(s_y, (h - surface(i, j + dj)) / dy);
    }
  }

  return std::sqrt(s_x * s_x + s_y * s_y);
}

} // end of namespace terminus
} // end of namespace debris
} // end of namespace pism
