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

// Credits
// -------
// ray_horizon() and surface_normal() are C++ re-implementations of algorithms from the
// "solshade" package by Aman Chokshi (https://github.com/amanchokshi/solshade, MIT License,
// (c) 2025 Aman Chokshi; Chokshi et al., Journal of Open Source Software,
// doi:10.21105/joss.09944): compute_horizon_map / _compute_horizon and
// compute_slope_aspect_normals (solshade/terrain.py).
//
// sky_view_factor() implements Dozier & Frew (1990)
//
// sun_position() is standard topocentric solar geometry; only the East-North-Up sun-vector
// convention follows solshade (solshade/solar.py).

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

// Adapted from solshade's compute_horizon_map / _compute_horizon (solshade/terrain.py).
// See the per-file credits at the top of this file.
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

// Adapted from solshade's compute_slope_aspect_normals (solshade/terrain.py): equivalent to
// its slope/aspect -> ENU normal, written directly in terms of the gradients. See the
// per-file credits at the top of this file.
void surface_normal(double dzdE, double dzdN, double &nE, double &nN, double &nU) {
  // Upward normal to the surface z = f(E, N) is proportional to (-f_E, -f_N, 1).
  double norm = std::sqrt(dzdE * dzdE + dzdN * dzdN + 1.0);
  nE = -dzdE / norm;
  nN = -dzdN / norm;
  nU = 1.0 / norm;
}

// Standard topocentric solar geometry (textbook spherical astronomy). The ENU sun-vector
// convention (E = cos(alt) sin(az), N = cos(alt) cos(az), U = sin(alt)) matches solshade's
// solar.py; the altitude/azimuth formulas themselves are standard.
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

// Implements the slope-corrected sky-view factor of Dozier & Frew (1990), "Rapid
// calculation of terrain parameters for radiation modeling from digital elevation data",
// IEEE Trans. Geosci. Remote Sens. 28(5):963-969, doi:10.1109/36.58986. The same
// formulation is used by TopoCalc (https://github.com/USDA-ARS-NWRC/topocalc). Not from
// solshade.
double sky_view_factor(const double *horizon, const double *azimuth, int n_dir,
                       double slope, double aspect) {
  const double cs = std::cos(slope);
  const double ss = std::sin(slope);

  double acc = 0.0;
  for (int k = 0; k < n_dir; ++k) {
    // Dozier & Frew use the horizon measured from the zenith. A terrain horizon below the
    // horizontal (negative elevation, e.g. on a peak) adds no sky, so clamp at the
    // horizontal.
    double h = horizon[k] > 0.0 ? horizon[k] : 0.0;
    double Hz = 0.5 * M_PI - h; // zenith angle of the visible-sky edge

    acc += cs * std::sin(Hz) * std::sin(Hz) +
           ss * std::cos(azimuth[k] - aspect) * (Hz - std::sin(Hz) * std::cos(Hz));
  }

  double svf = acc / n_dir; // (1/2pi) * integral, with d(azimuth) = 2pi/n_dir
  return clip(svf, 0.0, 1.0);
}

} // end of namespace terrain
} // end of namespace surface
} // end of namespace pism
