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

/*!
 * Unit tests for the terrain-horizon and sky-view kernels used by the "dEBM-enhanced"
 * surface model (see `pism/coupler/surface/terrain_insolation_kernel.hh`).
 *
 * These kernels are free of PISM, PETSc, and MPI, so this test links nothing but the
 * kernel itself and can be run directly:
 *
 *     ./pism_test_terrain_insolation_kernel
 *
 * Every case below is an analytic result: a linear DEM is reproduced exactly by bilinear
 * sampling, a constant slope has a constant horizon angle, solar noon puts the Sun due
 * south at a known altitude, and a flat site under a uniform horizon `h` has a sky-view
 * factor of cos(h)^2.
 */

#include "pism/coupler/surface/terrain_insolation_kernel.hh"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace pism::surface::terrain;

namespace {

const double pi = 3.14159265358979323846;

int failures = 0;

void check_close(const char *what, double actual, double expected, double tol) {
  if (not(std::fabs(actual - expected) <= tol)) {
    std::fprintf(stderr, "FAILED: %s: got %.15g, expected %.15g (tolerance %g)\n", what,
                 actual, expected, tol);
    ++failures;
  }
}

void check_true(const char *what, bool ok) {
  if (not ok) {
    std::fprintf(stderr, "FAILED: %s\n", what);
    ++failures;
  }
}

// A DEM tilted by `slope` (rise over run) toward the east: z = slope * x. Bilinear
// sampling reproduces it exactly, so the horizon angles below are exact.
std::vector<double> tilted_dem(int Mx, int My, double dx, double slope) {
  std::vector<double> dem(Mx * My);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      dem[j * Mx + i] = slope * (i * dx);
    }
  }
  return dem;
}

void test_sample_bilinear() {
  const int Mx = 5, My = 4;
  // z = 2 * i + 10 * j: linear, so bilinear interpolation is exact everywhere
  std::vector<double> dem(Mx * My);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      dem[j * Mx + i] = 2.0 * i + 10.0 * j;
    }
  }

  // exact at grid nodes
  check_close("bilinear at (0, 0)", sample_bilinear(dem.data(), Mx, My, 0, 0), 0.0, 1e-12);
  check_close("bilinear at (4, 3)", sample_bilinear(dem.data(), Mx, My, 4, 3), 38.0, 1e-12);

  // exact between nodes, because the DEM is linear
  check_close("bilinear at (1.5, 0)", sample_bilinear(dem.data(), Mx, My, 1.5, 0), 3.0, 1e-12);
  check_close("bilinear at (0, 2.25)", sample_bilinear(dem.data(), Mx, My, 0, 2.25), 22.5, 1e-12);
  check_close("bilinear at (2.5, 1.5)", sample_bilinear(dem.data(), Mx, My, 2.5, 1.5), 20.0, 1e-12);

  // coordinates outside the grid are clamped to the edge
  check_close("bilinear clamps fi < 0", sample_bilinear(dem.data(), Mx, My, -3.0, 1),
              sample_bilinear(dem.data(), Mx, My, 0.0, 1), 1e-12);
  check_close("bilinear clamps fi > Mx-1", sample_bilinear(dem.data(), Mx, My, 7.0, 1),
              sample_bilinear(dem.data(), Mx, My, 4.0, 1), 1e-12);
  check_close("bilinear clamps fj < 0", sample_bilinear(dem.data(), Mx, My, 2, -1.0),
              sample_bilinear(dem.data(), Mx, My, 2, 0.0), 1e-12);
  check_close("bilinear clamps fj > My-1", sample_bilinear(dem.data(), Mx, My, 2, 9.0),
              sample_bilinear(dem.data(), Mx, My, 2, 3.0), 1e-12);
}

void test_ray_horizon_flat() {
  const int Mx = 21, My = 21;
  const double dx = 100.0, dy = 100.0;
  std::vector<double> dem(Mx * My, 1234.5); // flat, but not at zero

  // A flat DEM has a horizon of exactly zero in every direction.
  for (int k = 0; k < 8; ++k) {
    double azimuth = 2.0 * pi * k / 8.0;
    check_close("flat DEM: horizon is zero",
                ray_horizon(dem.data(), Mx, My, dx, dy, 10, 10, azimuth, 50.0, 500.0), 0.0,
                1e-12);
  }
}

void test_ray_horizon_constant_slope() {
  const int Mx = 41, My = 41;
  const double dx = 100.0, dy = 100.0;
  const double slope = 0.1;
  auto dem = tilted_dem(Mx, My, dx, slope);

  const int i0 = 20, j0 = 20;
  const double step = 50.0, max_distance = 1000.0;

  // Looking upslope (east): the ray sees a constant angle atan(slope) at every distance.
  check_close("constant slope: horizon looking east",
              ray_horizon(dem.data(), Mx, My, dx, dy, i0, j0, 0.5 * pi, step, max_distance),
              std::atan(slope), 1e-12);

  // Looking downslope (west): the terrain falls away, so the "horizon" is negative.
  check_close("constant slope: horizon looking west",
              ray_horizon(dem.data(), Mx, My, dx, dy, i0, j0, 1.5 * pi, step, max_distance),
              -std::atan(slope), 1e-12);

  // Along strike (north and south) the elevation does not change.
  check_close("constant slope: horizon looking north",
              ray_horizon(dem.data(), Mx, My, dx, dy, i0, j0, 0.0, step, max_distance), 0.0,
              1e-12);
  check_close("constant slope: horizon looking south",
              ray_horizon(dem.data(), Mx, My, dx, dy, i0, j0, pi, step, max_distance), 0.0,
              1e-9);
}

void test_ray_horizon_wall() {
  const int Mx = 41, My = 41;
  const double dx = 100.0, dy = 100.0;

  // Flat ground with a 500 m wall filling every cell with i >= 25.
  std::vector<double> dem(Mx * My, 0.0);
  for (int j = 0; j < My; ++j) {
    for (int i = 25; i < Mx; ++i) {
      dem[j * Mx + i] = 500.0;
    }
  }

  // From i0 = 20 the foot of the wall is 5 cells = 500 m to the east, and it is 500 m
  // tall, so the horizon there is exactly atan2(500, 500) = pi/4. Points further along
  // the ray are the same height but further away, so they subtend a smaller angle.
  check_close("wall: horizon looking east",
              ray_horizon(dem.data(), Mx, My, dx, dy, 20, 20, 0.5 * pi, 50.0, 1500.0),
              0.25 * pi, 1e-12);

  // The wall is behind the observer when looking west.
  check_close("wall: horizon looking west",
              ray_horizon(dem.data(), Mx, My, dx, dy, 20, 20, 1.5 * pi, 50.0, 1500.0), 0.0,
              1e-12);
}

void test_ray_horizon_boundary() {
  const int Mx = 41, My = 41;
  const double dx = 100.0, dy = 100.0;
  auto dem = tilted_dem(Mx, My, dx, 0.1);

  // A ray leaving the domain on its very first step returns 0 even though the terrain
  // does slope away to the west (documented behavior: rays are clamped at the boundary).
  check_close("boundary cell: ray leaves the domain immediately",
              ray_horizon(dem.data(), Mx, My, dx, dy, 0, 20, 1.5 * pi, 50.0, 1000.0), 0.0,
              1e-12);
  check_close("boundary cell: ray leaves the domain immediately (south)",
              ray_horizon(dem.data(), Mx, My, dx, dy, 20, 0, pi, 50.0, 1000.0), 0.0, 1e-12);

  // A ray that exits through one boundary must not keep marching along the edge. Here the
  // ground is flat except for a tall block sitting on the northern edge, far to the east.
  // A northeast-bound ray from (20, 35) leaves the domain (fj > My - 1) after about 707 m,
  // by which point it has only reached i = 25. If the ray were merely clamped to the edge
  // instead of stopped, it would slide east along j = My - 1, run into the block at about
  // 2121 m, and report a horizon of roughly atan2(5000, 2121) = 1.17 radians.
  std::vector<double> edge(Mx * My, 0.0);
  for (int i = 33; i <= 37; ++i) {
    edge[(My - 1) * Mx + i] = 5000.0;
  }

  check_close("a ray that exits the domain does not slide along the boundary",
              ray_horizon(edge.data(), Mx, My, dx, dy, 20, 35, 0.25 * pi, 25.0, 3000.0), 0.0,
              1e-12);
}

void test_surface_normal() {
  double nE = 0.0, nN = 0.0, nU = 0.0;

  // A horizontal surface has a straight-up normal.
  surface_normal(0.0, 0.0, nE, nN, nU);
  check_close("flat surface normal, East", nE, 0.0, 1e-15);
  check_close("flat surface normal, North", nN, 0.0, 1e-15);
  check_close("flat surface normal, Up", nU, 1.0, 1e-15);

  // A 45-degree slope rising toward the east: the upward normal leans west.
  surface_normal(1.0, 0.0, nE, nN, nU);
  check_close("45-degree east slope, East", nE, -1.0 / std::sqrt(2.0), 1e-15);
  check_close("45-degree east slope, North", nN, 0.0, 1e-15);
  check_close("45-degree east slope, Up", nU, 1.0 / std::sqrt(2.0), 1e-15);

  // The normal is always a unit vector pointing up, and it leans opposite the gradient.
  const double gradients[5] = { -2.0, -0.5, 0.0, 0.25, 3.0 };
  for (double dzdE : gradients) {
    for (double dzdN : gradients) {
      surface_normal(dzdE, dzdN, nE, nN, nU);

      check_close("surface normal is a unit vector",
                  std::sqrt(nE * nE + nN * nN + nU * nU), 1.0, 1e-14);
      check_true("surface normal points up", nU > 0.0);
      check_true("surface normal leans away from the gradient",
                 nE * dzdE <= 0.0 and nN * dzdN <= 0.0);
    }
  }
}

void test_sun_position() {
  double altitude = 0.0, azimuth = 0.0;

  // Mid-northern latitude at equinox, solar noon: the Sun is due south, at an altitude of
  // 90 degrees minus the latitude.
  sun_position(0.25 * pi, 0.0, 0.0, altitude, azimuth);
  check_close("equinox noon at 45N: altitude", altitude, 0.25 * pi, 1e-12);
  check_close("equinox noon at 45N: azimuth is due south", azimuth, pi, 1e-12);

  // At the equator on the equinox the Sun rises due east and sets due west.
  sun_position(0.0, 0.0, -0.5 * pi, altitude, azimuth);
  check_close("equinox sunrise at the equator: altitude", altitude, 0.0, 1e-12);
  check_close("equinox sunrise at the equator: azimuth is due east", azimuth, 0.5 * pi, 1e-12);

  sun_position(0.0, 0.0, 0.5 * pi, altitude, azimuth);
  check_close("equinox sunset at the equator: altitude", altitude, 0.0, 1e-12);
  check_close("equinox sunset at the equator: azimuth is due west", azimuth, 1.5 * pi, 1e-12);

  // At solar noon the altitude is pi/2 - |latitude - declination|.
  const double latitudes[4]    = { -1.2, -0.3, 0.4, 1.1 };
  const double declinations[3] = { -0.409, 0.0, 0.409 };
  for (double latitude : latitudes) {
    for (double declination : declinations) {
      sun_position(latitude, declination, 0.0, altitude, azimuth);
      check_close("solar noon altitude", altitude,
                  0.5 * pi - std::fabs(latitude - declination), 1e-12);
    }
  }

  // Degenerate geometry returns a zero azimuth rather than a NaN: the Sun at the zenith,
  // and an observer at the north pole.
  sun_position(0.0, 0.0, 0.0, altitude, azimuth);
  check_close("Sun at the zenith: altitude", altitude, 0.5 * pi, 1e-12);
  check_close("Sun at the zenith: azimuth", azimuth, 0.0, 1e-15);

  sun_position(0.5 * pi, 0.2, 1.0, altitude, azimuth);
  check_true("north pole: altitude is finite", std::isfinite(altitude));
  check_close("north pole: azimuth", azimuth, 0.0, 1e-15);

  // Over a sweep of geometries the outputs stay in range, and the altitude is symmetric
  // about solar noon (morning and afternoon mirror each other).
  for (int a = 0; a < 17; ++a) {
    double latitude = -0.5 * pi + pi * a / 16.0;
    for (int b = 0; b < 9; ++b) {
      double declination = -0.409 + 0.818 * b / 8.0;
      for (int c = 0; c < 25; ++c) {
        double hour_angle = -pi + 2.0 * pi * c / 24.0;

        sun_position(latitude, declination, hour_angle, altitude, azimuth);

        check_true("altitude is in [-pi/2, pi/2]",
                   altitude >= -0.5 * pi - 1e-12 and altitude <= 0.5 * pi + 1e-12);
        check_true("azimuth is in [0, 2pi)", azimuth >= 0.0 and azimuth < 2.0 * pi);

        double mirrored = 0.0, unused = 0.0;
        sun_position(latitude, declination, -hour_angle, mirrored, unused);
        check_close("altitude is symmetric about solar noon", altitude, mirrored, 1e-12);
      }
    }
  }
}

// Equally spaced azimuths, clockwise from north, as used by TerrainInsolation.
std::vector<double> azimuths(int n_dir) {
  std::vector<double> result(n_dir);
  for (int k = 0; k < n_dir; ++k) {
    result[k] = 2.0 * pi * k / n_dir;
  }
  return result;
}

void test_sky_view_factor() {
  const int n_dir = 72;
  auto azimuth = azimuths(n_dir);

  // A horizontal site with nothing blocking the sky sees all of it.
  {
    std::vector<double> horizon(n_dir, 0.0);
    check_close("horizontal and unobstructed: sky view is 1",
                sky_view_factor(horizon.data(), azimuth.data(), n_dir, 0.0, 0.0), 1.0, 1e-12);
  }

  // A horizontal site ringed by a uniform horizon `h` sees cos(h)^2 of the sky.
  const double heights[5] = { 0.1, 0.3, 0.5, 0.9, 1.4 };
  for (double h : heights) {
    std::vector<double> horizon(n_dir, h);
    check_close("horizontal under a uniform horizon: sky view is cos(h)^2",
                sky_view_factor(horizon.data(), azimuth.data(), n_dir, 0.0, 0.0),
                std::cos(h) * std::cos(h), 1e-12);
  }

  // A horizon *below* the horizontal (a site on a peak) adds no sky beyond the hemisphere.
  {
    std::vector<double> horizon(n_dir, -0.5);
    check_close("a horizon below the horizontal is clamped",
                sky_view_factor(horizon.data(), azimuth.data(), n_dir, 0.0, 0.0), 1.0, 1e-12);
  }

  // Raising a uniform horizon can only take sky away.
  {
    double previous = 2.0;
    for (int k = 0; k <= 15; ++k) {
      std::vector<double> horizon(n_dir, 0.5 * pi * k / 15.0);
      double svf = sky_view_factor(horizon.data(), azimuth.data(), n_dir, 0.0, 0.0);
      check_true("sky view decreases as a uniform horizon rises", svf < previous + 1e-12);
      check_true("sky view is in [0, 1]", svf >= 0.0 and svf <= 1.0);
      previous = svf;
    }
  }

  // Tilting an unobstructed surface: the Dozier & Frew slope term averages to zero over a
  // full circle of azimuths, leaving cos(slope).
  for (double slope : { 0.0, 0.2, 0.5, 1.0 }) {
    std::vector<double> horizon(n_dir, 0.0);
    check_close("tilted and unobstructed: sky view is cos(slope)",
                sky_view_factor(horizon.data(), azimuth.data(), n_dir, slope, 0.7),
                std::cos(slope), 1e-12);
  }

  // The slope term is direction-dependent: an obstruction downslope (in the direction the
  // surface faces) and the same obstruction upslope do not have the same effect. This
  // pins down the sin(slope) * cos(azimuth - aspect) term, which the tests above cannot
  // see because it integrates away.
  {
    const double slope = 0.5, aspect = 0.0, h = 1.0;

    std::vector<double> downslope(n_dir, 0.0), upslope(n_dir, 0.0);
    for (int k = 0; k < n_dir; ++k) {
      // block a quadrant centered on the aspect, and the opposite quadrant
      double da = std::fabs(std::remainder(azimuth[k] - aspect, 2.0 * pi));
      if (da < 0.25 * pi) {
        downslope[k] = h;
      }
      if (da > 0.75 * pi) {
        upslope[k] = h;
      }
    }

    double svf_down = sky_view_factor(downslope.data(), azimuth.data(), n_dir, slope, aspect);
    double svf_up   = sky_view_factor(upslope.data(), azimuth.data(), n_dir, slope, aspect);

    check_true("obstruction up- and downslope differ", std::fabs(svf_down - svf_up) > 0.1);
    check_true("blocking the downslope sky removes more than blocking upslope",
               svf_down < svf_up);
    check_true("both are in [0, 1]",
               svf_down >= 0.0 and svf_down <= 1.0 and svf_up >= 0.0 and svf_up <= 1.0);
  }

  // The result stays in [0, 1] even for geometries that push the formula out of range.
  for (int a = 0; a <= 10; ++a) {
    double slope = 0.5 * pi * a / 10.0;
    for (int b = 0; b <= 10; ++b) {
      std::vector<double> horizon(n_dir, 0.5 * pi * b / 10.0);
      double svf = sky_view_factor(horizon.data(), azimuth.data(), n_dir, slope, 1.3);
      check_true("sky view stays in [0, 1]", svf >= 0.0 and svf <= 1.0);
    }
  }
}

} // end of anonymous namespace

int main() {
  test_sample_bilinear();
  test_ray_horizon_flat();
  test_ray_horizon_constant_slope();
  test_ray_horizon_wall();
  test_ray_horizon_boundary();
  test_surface_normal();
  test_sun_position();
  test_sky_view_factor();

  if (failures > 0) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }

  std::printf("terrain insolation kernel: all checks passed\n");
  return 0;
}
