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

#include <cmath>
#include <limits>

namespace pism {
namespace surface {
namespace terrain {

static inline double clip(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

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

double ray_horizon(const double *dem, int Mx, int My, double dx, double dy,
                   int i0, int j0, double azimuth, double step, double max_distance) {
  const double z0 = dem[j0 * Mx + i0];

  // Azimuth clockwise from north: horizontal direction (East, North) = (sin, cos).
  const double sE = std::sin(azimuth);
  const double cN = std::cos(azimuth);

  double best = -std::numeric_limits<double>::infinity();

  for (double d = step; d <= max_distance; d += step) {
    double fi = i0 + sE * d / dx;
    double fj = j0 + cN * d / dy;

    // Stop the ray once it leaves the (non-periodic) physical domain.
    if (fi < 0.0 || fi > Mx - 1 || fj < 0.0 || fj > My - 1) {
      break;
    }

    double z = sample_bilinear(dem, Mx, My, fi, fj);
    double angle = std::atan2(z - z0, d);
    if (angle > best) {
      best = angle;
    }
  }

  return std::isfinite(best) ? best : 0.0;
}

void surface_normal(double dzdE, double dzdN, double &nE, double &nN, double &nU) {
  // Upward normal to the surface z = f(E, N) is proportional to (-f_E, -f_N, 1).
  double norm = std::sqrt(dzdE * dzdE + dzdN * dzdN + 1.0);
  nE = -dzdE / norm;
  nN = -dzdN / norm;
  nU = 1.0 / norm;
}

void sun_position(double latitude, double declination, double hour_angle,
                  double &altitude, double &azimuth) {
  const double sl = std::sin(latitude), cl = std::cos(latitude);
  const double sd = std::sin(declination), cd = std::cos(declination);
  const double sH = std::sin(hour_angle), cH = std::cos(hour_angle);

  double sin_alt = clip(sl * sd + cl * cd * cH, -1.0, 1.0);
  altitude = std::asin(sin_alt);

  double cos_alt = std::cos(altitude);

  // Degenerate geometry: sun at the zenith, or observer at a geographic pole. Azimuth is
  // undefined; return 0 (irrelevant for the cosine projection at the zenith, and PISM
  // domains are not located exactly at a pole).
  if (cos_alt < 1e-8 || cl < 1e-8) {
    azimuth = 0.0;
    return;
  }

  double sinA = -cd * sH / cos_alt;
  double cosA = (sd - sl * sin_alt) / (cl * cos_alt);

  double A = std::atan2(sinA, cosA); // clockwise from north, in (-pi, pi]
  if (A < 0.0) {
    A += 2.0 * M_PI;
  }
  azimuth = A;
}

} // end of namespace terrain
} // end of namespace surface
} // end of namespace pism
