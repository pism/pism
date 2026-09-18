// Copyright (C) 2026 PISM Authors
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

#include "pism/coupler/surface/terrain_insolation_kernel.hh"
#include "pism/util/SunPosition.hh"

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

// References
// ----------
//
// ray_horizon() is a C++ re-implementations of the algorithm from the "solshade" package
// by Aman Chokshi (https://github.com/amanchokshi/solshade, MIT License, (c) 2025 Aman
// Chokshi;
//
// Chokshi et al., Journal of Open Source Software, doi:10.21105/joss.09944):
// compute_horizon_map / _compute_horizon and compute_slope_aspect_normals
// (solshade/terrain.py).
//
// J. Dozier and J. Frew, “Rapid calculation of terrain parameters for radiation modeling
// from digital elevation data,” IEEE Transactions on Geoscience and Remote Sensing, vol.
// 28, no. 5, pp. 963–969, 1990, doi: 10.1109/36.58986.

namespace pism {
namespace surface {
namespace terrain {

static inline double clip(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

/*!
 * Sample a row-major array `dem` with size `Mx` columns by `My` rows at location
 * `fi`,`fj` ("fractional" indexes in X and Y directions) using bilinear interpolation.
 */
double sample_bilinear(const double *dem, int Mx, int My, double fi, double fj) {
  // Clamp to the valid index range (samples on or just outside the boundary fall back to
  // the edge value). This mirrors the clamping in Grid::compute_point_neighbors.
  fi = clip(fi, 0.0, (double)(Mx - 1));
  fj = clip(fj, 0.0, (double)(My - 1));

  int i0 = (int)std::floor(fi);
  int j0 = (int)std::floor(fj);
  int i1 = i0 < Mx - 1 ? i0 + 1 : i0;
  int j1 = j0 < My - 1 ? j0 + 1 : j0;

  double a = fi - i0; // weight toward i1
  double b = fj - j0; // weight toward j1

  double z00 = dem[j0 * Mx + i0];
  double z10 = dem[j0 * Mx + i1];
  double z01 = dem[j1 * Mx + i0];
  double z11 = dem[j1 * Mx + i1];

  return (1.0 - a) * (1.0 - b) * z00 + a * (1.0 - b) * z10 +
         (1.0 - a) * b * z01 + a * b * z11;
}

/*!
 * Compute the horizon map at a grid point.
 *
 * @param[in] dem row-major array containing the surface elevation DEM
 * @param[in] Mx grid size in the X direction
 * @param[in] My grid size in the Y direction
 * @param[in] dx grid spacing in the X direction
 * @param[in] dy grid spacing in the Y direction
 * @param[in] i0 X-index of a grid point
 * @param[in] j0 Y-index of a grid point
 * @param[in] azimuth direction (clockwise from the direction of the Y axis) to consider
 * @param[in] step step length (meters)
 * @param[in] max_distance maximum distance to consider (meters)
 *
 * Returns the altitude angle of the horizon in direction corresponding to `azimuth`, in
 * radian.
 *
 * Adapted from solshade's compute_horizon_map / _compute_horizon (solshade/terrain.py).
 */
double ray_horizon(const double *dem, int Mx, int My, double dx, double dy,
                   int i0, int j0, double azimuth, double step, double max_distance) {
  const double z0 = dem[j0 * Mx + i0];

  assert(step > 0.0);

  // Azimuth A is measured clockwise from the Y direction, so the (x,y) vector in the
  // direction A has the form (x, y) = (sin(A), cos(A)).
  const double v_x = std::sin(azimuth);
  const double v_y = std::cos(azimuth);

  double max_slope = -std::numeric_limits<double>::infinity();

  for (double d = step; d <= max_distance; d += step) {
    double fi = i0 + v_x * d / dx;
    double fj = j0 + v_y * d / dy;

    // Stop the ray once it leaves the (non-periodic) physical domain.
    if (fi < 0.0 || fi > Mx - 1 || fj < 0.0 || fj > My - 1) {
      break;
    }

    double z = sample_bilinear(dem, Mx, My, fi, fj);

    max_slope = std::fmax(max_slope, (z - z0) / d);
  }

  return std::isfinite(max_slope) ? std::atan(max_slope) : 0.0;
}

void sun_position(double latitude, double declination, double hour_angle,
                  double &altitude, double &azimuth) {
  SunPosition sp(declination);
  sp.set_latitude(latitude);
  sp.compute(hour_angle, altitude, azimuth);
}

/*!
 * Implements the slope-corrected sky-view factor of Dozier & Frew (1990).
 *
 * @param[in] horizon array of `n_dir` elements containing horizon altitudes in directions
 *                    in `azimuth`, in radians
 * @param[in] azimuth array of `n_dir` elements containing azimuth directions, in radians
 * @param[in] n_dir lengths of arrays `horizon` and `azimuth`
 * @param[in] slope slope at the current location
 * @param[in] aspect aspect at the current location (radians, clockwise from the north)
 *
 * Returns the sky view factor between 0 and 1.
 */
double sky_view_factor(const double *horizon, const double *azimuth, int n_dir,
                       double slope, double aspect) {

  double cos_slope = std::cos(slope);
  double sin_slope = std::sin(slope);

  // Note: this code implements equation (7b) in Dozier and Frew.
  double acc = 0.0;
  for (int k = 0; k < n_dir; ++k) {
    // Dozier & Frew use the horizon measured from the zenith. A terrain horizon below the
    // horizontal (negative elevation, e.g. at a peak) adds no sky, so clamp at the
    // horizontal.
    double Hz = M_PI_2 - std::fmax(horizon[k], 0.0); // zenith angle of the visible-sky edge
    double sin_Hz = std::sin(Hz);

    acc += cos_slope * sin_Hz * sin_Hz +
           sin_slope * std::cos(azimuth[k] - aspect) * (Hz - sin_Hz * std::cos(Hz));
  }

  double svf = acc / n_dir; // (1/2pi) * integral, with d(azimuth) = 2pi/n_dir
  return clip(svf, 0.0, 1.0);
}

} // end of namespace terrain
} // end of namespace surface
} // end of namespace pism
