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

#include "pism/coupler/surface/TerrainInsolation.hh"
#include "pism/coupler/surface/terrain_insolation_kernel.hh"

#include <cmath>
#include <geodesic.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <vector>

#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/Logger.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/LonLatCalculator.hh"
#include "pism/util/SunPosition.hh"

#if (Pism_USE_PROJ == 0)
#error "This code requires PROJ"
#endif
#include <proj.h>

/*!
 * Reference:
 *
 * A. Chokshi, “Solshade: Terrain-aware Solar Illumination Modelling using Digital
 * Elevation Models and Orbital Geometry,” Journal of Open Source Software, vol. 11, no.
 * 122, p. 9944, Jun. 2026, doi: 10.21105/joss.09944.
 */

namespace pism {
namespace surface {

TerrainInsolation::TerrainInsolation(std::shared_ptr<const Grid> grid,
                                     std::function<double(double)> atmosphere_transmissivity)
    : m_grid(grid), m_insolation(grid, "insolation"),
      m_orbital_parameters(*grid->ctx()),
      m_transmissivity(atmosphere_transmissivity) {

  m_insolation.metadata(0)
      .long_name("daily mean terrain-shaded surface insolation")
      .units("W m^-2");

  auto config = m_grid->ctx()->config();

  m_n_directions  = static_cast<int>(config->get_number("surface.debm_enhanced.horizon.n_directions"));
  m_max_distance  = config->get_number("surface.debm_enhanced.horizon.max_distance");
  m_step          = config->get_number("surface.debm_enhanced.horizon.step");
  m_insolation_dt  = config->get_number("surface.debm_enhanced.insolation_dt");
  m_solar_constant = config->get_number("surface.debm_simple.solar_constant");
  m_diffuse_fraction = config->get_number("surface.debm_enhanced.diffuse_fraction");

  if (not (m_step > 0.0)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "surface.debm_enhanced.horizon.step must be positive");
  }

  if (m_n_directions < 1) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "surface.debm_enhanced.horizon.n_directions must be positive");
  }

  // azimuth sample directions (radians), clockwise from north
  std::vector<double> azimuth(m_n_directions);
  for (int k = 0; k < m_n_directions; ++k) {
    azimuth[k] = 2.0 * M_PI * k / m_n_directions;
  }

  m_horizon = std::make_shared<array::Array3D>(m_grid, "horizon", array::WITHOUT_GHOSTS,
                                               azimuth);
  m_horizon->metadata(0)
      .long_name("terrain horizon elevation angle as a function of azimuth")
      .units("radian");

  bool use_sky_view  = config->get_flag("surface.debm_enhanced.use_sky_view_factor");
  if (use_sky_view) {
    m_sky_view = std::make_shared<array::Scalar>(m_grid, "sky_view_factor");
    m_sky_view->metadata(0)
        .long_name("sky-view factor (fraction of the diffuse sky hemisphere visible "
                   "from the terrain-shaded, tilted surface)")
        .units("1");
  }

  // Allocate the scatter to all ranks and the vector that will hold local copies of the
  // DEM:
  PetscErrorCode ierr;
  ierr = DMDAGlobalToNaturalAllCreate(*m_grid->get_dm(1, 0), m_scatter.rawptr());
  PISM_CHK(ierr, "DMDAGlobalToNaturalAllCreate");

  ierr = VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(m_grid->Mx() * m_grid->My()),
                      m_dem_local.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");
}

const array::Scalar& TerrainInsolation::insolation() const {
  return m_insolation;
}

const array::Array3D &TerrainInsolation::horizon() const {
  return *m_horizon;
}

const array::Scalar &TerrainInsolation::sky_view() const {
  if (not sky_view_enabled()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "sky view factor is not available");
  }
  return *m_sky_view;
}

bool TerrainInsolation::sky_view_enabled() const {
  return m_sky_view != nullptr;
}

/*!
 * Compute the angle between the Y axis of the grid (vector (0, 1)) and the meridian
 * (vector pointing from the current point to the north pole), in radians.
 *
 * See https://proj.org/en/stable/geodesic.html
 */
void compute_azimuth(array::Scalar &output) {
  auto grid = output.grid();

  std::string proj_string = grid->get_mapping_info()["proj_params"];

  double a;
  double invf;
  // get ellipsoid information from PROJ:
  {
    PJ_CONTEXT *C = proj_context_create();
    PJ *P         = proj_create(C, proj_string.c_str());
    proj_ellipsoid_get_parameters(C, proj_get_ellipsoid(C, P), &a, 0, 0, &invf);
    proj_destroy(P);
    proj_context_destroy(C);
  }

  LonLatCalculator calc(proj_string);

  geod_geodesic g;
  geod_init(&g, a, invf != 0 ? 1 / invf : 0);

  double dx = grid->dx();

  array::AccessScope list{&output};
  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    double x = grid->x(i);
    double y = grid->y(j);

    auto pt = calc.lonlat(x, y);
    double lon1 = pt[0], lat1 = pt[1];

    pt = calc.lonlat(x + dx, y);
    double lon2 = pt[0], lat2 = pt[1];

    double A = 0.0;
    geod_inverse(&g, lat1, lon1, lat2, lon2, nullptr, &A, nullptr);

    output(i, j) = proj_torad(90.0 - A);
  }
}

void TerrainInsolation::update_horizon_map(const array::Scalar1 &surface_elevation) {
  auto log = m_grid->ctx()->log();

  log->message(2, "* Updating the horizon map...\n");
  double start = get_time(m_grid->com);

  const auto &profiling = m_grid->ctx()->profiling();

  const int Mx = static_cast<int>(m_grid->Mx());
  const int My = static_cast<int>(m_grid->My());
  const double dx = m_grid->dx();
  const double dy = m_grid->dy();

  auto diff_x = [Mx, dx](const array::Scalar1 &F, int i, int j) {
    // use one-sided finite differences at domain boundaries:
    int ip = i < Mx - 1 ? i + 1 : i;
    int im = i > 0 ? i - 1 : i;

    return (F(ip, j) - F(im, j)) / ((ip - im) * dx);
  };

  auto diff_y = [My, dy](const array::Scalar1 &F, int i, int j) {
    // use one-sided finite differences at domain boundaries:
    int jp = j < My - 1 ? j + 1 : j;
    int jm = j > 0 ? j - 1 : j;

    return (F(i, jp) - F(i, jm)) / ((jp - jm) * dy);
  };

  // Scatter the full DEM to every rank. After this block each rank has a copy of the
  // global surface elevation that lets it compute shading at its owned cells without
  // ghost communication.
  {
    profiling.begin("surface.debm_enhanced.scatter_dem");
    PetscErrorCode ierr;
    auto dm = surface_elevation.dm();
    petsc::TemporaryGlobalVec dem_global(dm);
    // Note: we use DMLocalToGlobal because surface_elevation is ghosted (local)
    ierr = DMLocalToGlobal(*dm, surface_elevation.vec(), INSERT_VALUES, dem_global);
    PISM_CHK(ierr, "DMLocalToGlobal");

    ierr = VecScatterBegin(m_scatter, dem_global, m_dem_local, INSERT_VALUES, SCATTER_FORWARD);
    PISM_CHK(ierr, "VecScatterBegin");
    ierr = VecScatterEnd(m_scatter, dem_global, m_dem_local, INSERT_VALUES, SCATTER_FORWARD);
    PISM_CHK(ierr, "VecScatterEnd");
    profiling.end("surface.debm_enhanced.scatter_dem");
  }

  // Compute surface normals (centered differences on the global DEM, one-sided at the
  // domain boundary) and the horizon map (the dominant cost) for every owned cell.
  profiling.begin("surface.debm_enhanced.horizon");

  const auto &azimuth = m_horizon->levels();

  petsc::VecArray dem(m_dem_local);
  array::AccessScope scope{ &surface_elevation, m_horizon.get() };

  bool use_sky_view = sky_view_enabled();
  if (use_sky_view) {
    scope.add(*m_sky_view);
  }

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    double *column = m_horizon->get_column(i, j);

    // Compute the upward-pointing normal to the surface:
    double s_x = diff_x(surface_elevation, i, j);
    double s_y = diff_y(surface_elevation, i, j);

    // FIXME: incorrect assumption!
    double s_e = s_x;
    double s_n = s_y;

    double Ne = -s_e;
    double Nn = -s_n;
    double Nu = 1.0;

    // Scale to get the unit normal:
    {
      // Note that norm != 0.0 because Nu == 1
      double norm = std::sqrt(Ne * Ne + Nn * Nn + Nu * Nu);

      Ne /= norm;
      Nn /= norm;
      Nu /= norm;
    }

    for (int k = 0; k < m_n_directions; ++k) {
      column[k] = terrain::ray_horizon(dem.get(), Mx, My, dx, dy, i, j, azimuth[k], m_step,
                                       m_max_distance);
    }

    // sky-view factor from the horizon and the surface slope/aspect (the latter recovered
    // from the unit normal: slope = acos(Nu), aspect = atan2(Ne, Nn), clockwise from north)
    if (use_sky_view) {
      double slope = std::acos(pism::clip(Nu, -1.0, 1.0));
      double aspect = std::atan2(Ne, Nn);
      (*m_sky_view)(i, j) =
          terrain::sky_view_factor(column, azimuth.data(), m_n_directions, slope,
                                   aspect);
    }
  }
  profiling.end("surface.debm_enhanced.horizon");
  double end = get_time(m_grid->com);
  log->message(2, "* Updated the horizon map in %f s.\n", end - start);
}

//! Periodic linear interpolation of a horizon column at the given azimuth (radians).
double TerrainInsolation::interpolate(const double *column, int n, double azimuth) {
  const double two_pi = 2.0 * M_PI;
  const double da = two_pi / n;

  double a = azimuth - two_pi * std::floor(azimuth / two_pi); // wrap to [0, 2*pi)
  double x = a / da;
  int k = static_cast<int>(std::floor(x));
  double frac = x - k;

  int k0 = k % n;
  int k1 = (k + 1) % n;

  return column[k0] * (1.0 - frac) + column[k1] * frac;
}

// The per-timestep shadow test, Lambertian cosine projection (max(0, normal . sun)), and
// inverse-square distance scaling are adapted from solshade's compute_flux_timeseries
// (solshade/irradiance.py).
//
// Here they are integrated over the diurnal cycle to compute daily energy.
void TerrainInsolation::update_daily_insolation(double time,
                                                const array::Scalar1 &surface_elevation) {

  auto p = m_orbital_parameters.compute(time);
  double declination = p.solar_declination;
  double distance_factor = p.distance_factor;

  const double seconds_per_day = 86400.0;

  SunPosition sun_position(declination);

  // number of sub-daily samples used to integrate the diurnal cycle
  int M = static_cast<int>(std::lround(seconds_per_day / m_insolation_dt));
  if (M < 1) {
    M = 1;
  }

  // Pre-compute cos() and sin() of M hour angles:
  {
    std::vector<double> hour_angles(M);
    for (int m = 0; m < M; ++m) {
      // hour angle sweeps the full day, midpoint rule over [-pi, pi)
      hour_angles[m] = -M_PI + 2.0 * M_PI * (m + 0.5) / M;
    }
    sun_position.set_hour_angles(hour_angles);
  }

  const double dt = seconds_per_day / M;

  const auto &profiling = m_grid->ctx()->profiling();
  profiling.begin("surface.debm_enhanced.daily_insolation");

  const auto &latitude = m_grid->latitude();

  array::AccessScope scope{ &latitude, &m_insolation, &surface_elevation, m_horizon.get() };

  bool use_sky_view = sky_view_enabled();
  if (use_sky_view) {
    scope.add(*m_sky_view);
  }

  int Mx = (int)m_grid->Mx();
  int My = (int)m_grid->My();
  double dx = m_grid->dx();
  double dy = m_grid->dy();

  auto diff_x = [Mx, dx](const array::Scalar1 &F, int i, int j) {
    // use one-sided finite differences at domain boundaries:
    int ip = i < Mx - 1 ? i + 1 : i;
    int im = i > 0 ? i - 1 : i;

    return (F(ip, j) - F(im, j)) / ((ip - im) * dx);
  };

  auto diff_y = [My, dy](const array::Scalar1 &F, int i, int j) {
    // use one-sided finite differences at domain boundaries:
    int jp = j < My - 1 ? j + 1 : j;
    int jm = j > 0 ? j - 1 : j;

    return (F(i, jp) - F(i, jm)) / ((jp - jm) * dy);
  };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    sun_position.set_latitude(latitude(i, j) * (M_PI / 180.0));

    // Compute the upward-pointing normal to the surface:
    double s_x = diff_x(surface_elevation, i, j);
    double s_y = diff_y(surface_elevation, i, j);

    // FIXME: incorrect assumption!
    double s_e = s_x;
    double s_n = s_y;

    double Ne = -s_e;
    double Nn = -s_n;
    double Nu = 1.0;

    // Scale to get the unit normal:
    {
      // Note that norm != 0.0 because Nu == 1
      double norm = std::sqrt(Ne * Ne + Nn * Nn + Nu * Nu);

      Ne /= norm;
      Nn /= norm;
      Nu /= norm;
    }

    const double *horizon = m_horizon->get_column(i, j);

    // Split into a direct-beam fraction (terrain-shaded) and an isotropic diffuse fraction
    // (reduced by the sky-view factor). With the sky-view factor disabled the diffuse term
    // is dropped and the result is pure direct beam.
    const double f_diff = use_sky_view ? m_diffuse_fraction : 0.0;
    const double svf = use_sky_view ? (*m_sky_view)(i, j) : 0.0;

    double energy = 0.0;
    // loop over hour angles:
    for (int hour_angle_idx = 0; hour_angle_idx < M; ++hour_angle_idx) {

      double altitude = 0.0, azimuth = 0.0;

      // vector pointing toward the center of the sun in the east-north-up coordinate
      // system:
      double solar_vector[3] = {0.0, 0.0, 0.0};

      sun_position.compute_at_set_hour_angle(hour_angle_idx, altitude, azimuth, solar_vector);

      if (altitude <= 0.0) {
        // Sun below the astronomical horizon: no direct and no diffuse contribution
        continue;
      }

      // Note: from equation (13) in Sproul2007 sin(altitude) is equal to the Z ("up")
      // component of the solar vector.
      double sin_altitude = solar_vector[2];
      // top-of-atmosphere horizontal irradiance, the basis for the diffuse component
      double toa_horizontal = m_solar_constant * distance_factor * sin_altitude;

      // diffuse: isotropic sky scaled by the sky-view factor; reaches shadowed cells too
      energy += f_diff * toa_horizontal * svf * dt;

      double altitude_threshold = interpolate(horizon, m_n_directions, azimuth);

      // direct beam: only when the Sun clears the local horizon and lights the surface
      if (altitude > altitude_threshold) {
        double Se = solar_vector[0];
        double Sn = solar_vector[1];
        double Su = solar_vector[2];

        double mu = Ne * Se + Nn * Sn + Nu * Su;
        if (mu > 0.0) {
          energy += (1.0 - f_diff) * m_solar_constant * distance_factor * mu * dt;
        }
      }
    } // end of the loop over hour angles

    // store the daily-mean insolation rate (W m-2), matching dEBM-simple's "insolation"
    // diagnostic units (the melt code multiplies this rate by the sub-step length)
    m_insolation(i, j) = m_transmissivity(surface_elevation(i, j)) * energy / seconds_per_day;
  }

  profiling.end("surface.debm_enhanced.daily_insolation");
}

} // end of namespace surface
} // end of namespace pism
