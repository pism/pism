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

#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {
namespace surface {

TerrainInsolation::TerrainInsolation(std::shared_ptr<const Grid> grid)
  : m_grid(grid) {

  auto config = m_grid->ctx()->config();

  m_n_directions  = static_cast<int>(config->get_number("surface.debm_enhanced.horizon.n_directions"));
  m_max_distance  = config->get_number("surface.debm_enhanced.horizon.max_distance");
  m_step          = config->get_number("surface.debm_enhanced.horizon.step");
  m_ephemeris_dt  = config->get_number("surface.debm_enhanced.horizon.ephemeris_dt");
  m_solar_constant = config->get_number("surface.debm_simple.solar_constant");
  m_use_sky_view  = config->get_flag("surface.debm_enhanced.use_sky_view_factor");

  if (m_n_directions < 1) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "surface.debm_enhanced.horizon.n_directions must be positive");
  }

  // azimuth sample directions (radians), clockwise from north
  m_azimuth.resize(m_n_directions);
  for (int k = 0; k < m_n_directions; ++k) {
    m_azimuth[k] = 2.0 * M_PI * k / m_n_directions;
  }

  m_horizon = std::make_shared<array::Array3D>(m_grid, "horizon", array::WITHOUT_GHOSTS,
                                               m_azimuth);
  m_horizon->metadata(0)
      .long_name("terrain horizon elevation angle as a function of azimuth")
      .units("radian");

  m_normal_e = std::make_shared<array::Scalar>(m_grid, "surface_normal_e");
  m_normal_n = std::make_shared<array::Scalar>(m_grid, "surface_normal_n");
  m_normal_u = std::make_shared<array::Scalar>(m_grid, "surface_normal_u");

  if (m_use_sky_view) {
    m_sky_view = std::make_shared<array::Scalar>(m_grid, "sky_view_factor");
    m_sky_view->metadata(0)
        .long_name("sky-view factor (fraction of the diffuse sky hemisphere visible "
                   "from the terrain-shaded, tilted surface)")
        .units("1");
  }
}

const array::Array3D &TerrainInsolation::horizon() const {
  return *m_horizon;
}

const array::Scalar &TerrainInsolation::sky_view() const {
  return *m_sky_view;
}

bool TerrainInsolation::sky_view_enabled() const {
  return m_use_sky_view;
}

void TerrainInsolation::init(const array::Scalar &surface_elevation) {
  const int Mx = static_cast<int>(m_grid->Mx());
  const int My = static_cast<int>(m_grid->My());
  const double dx = m_grid->dx();
  const double dy = m_grid->dy();

  // Gather the full DEM onto rank 0, then broadcast it to every rank. The result is a
  // contiguous, row-major (dem[j * Mx + i]) copy of the global surface elevation that lets
  // each rank ray-march its owned cells without any ghost communication.
  m_dem.resize(static_cast<size_t>(Mx) * static_cast<size_t>(My));

  auto work0 = surface_elevation.allocate_proc0_copy();
  surface_elevation.put_on_proc0(*work0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      petsc::VecArray array(*work0);
      const double *a = array.get();
      for (size_t k = 0; k < m_dem.size(); ++k) {
        m_dem[k] = a[k];
      }
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  int ierr = MPI_Bcast(m_dem.data(), static_cast<int>(m_dem.size()), MPI_DOUBLE, 0,
                       m_grid->com);
  PISM_CHK(ierr, "MPI_Bcast");

  // Compute surface normals (centred differences on the global DEM, one-sided at the
  // domain boundary) and the horizon map for every owned cell.
  const double *dem = m_dem.data();

  std::vector<double> column(m_n_directions);

  array::AccessScope scope{ m_normal_e.get(), m_normal_n.get(), m_normal_u.get(),
                            m_horizon.get() };
  if (m_use_sky_view) {
    scope.add(*m_sky_view);
  }

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    int ip = i < Mx - 1 ? i + 1 : i;
    int im = i > 0 ? i - 1 : i;
    int jp = j < My - 1 ? j + 1 : j;
    int jm = j > 0 ? j - 1 : j;

    double dzdE = (dem[j * Mx + ip] - dem[j * Mx + im]) / ((ip - im) * dx);
    double dzdN = (dem[jp * Mx + i] - dem[jm * Mx + i]) / ((jp - jm) * dy);

    double nE = 0.0, nN = 0.0, nU = 1.0;
    terrain::surface_normal(dzdE, dzdN, nE, nN, nU);
    (*m_normal_e)(i, j) = nE;
    (*m_normal_n)(i, j) = nN;
    (*m_normal_u)(i, j) = nU;

    for (int k = 0; k < m_n_directions; ++k) {
      column[k] = terrain::ray_horizon(dem, Mx, My, dx, dy, i, j, m_azimuth[k], m_step,
                                       m_max_distance);
    }
    m_horizon->set_column(i, j, column.data());

    // sky-view factor from the horizon and the surface slope/aspect (the latter recovered
    // from the unit normal: slope = acos(nU), aspect = atan2(nE, nN), clockwise from north)
    if (m_use_sky_view) {
      double slope = std::acos(nU < -1.0 ? -1.0 : (nU > 1.0 ? 1.0 : nU));
      double aspect = std::atan2(nE, nN);
      (*m_sky_view)(i, j) =
          terrain::sky_view_factor(column.data(), m_azimuth.data(), m_n_directions, slope,
                                   aspect);
    }
  }
}

//! Periodic linear interpolation of a horizon column at the given azimuth (radians).
double TerrainInsolation::horizon_at(const double *column, double azimuth) const {
  const double two_pi = 2.0 * M_PI;
  const double da = two_pi / m_n_directions;

  double a = azimuth - two_pi * std::floor(azimuth / two_pi); // wrap to [0, 2*pi)
  double x = a / da;
  int k = static_cast<int>(std::floor(x));
  double frac = x - k;

  int k0 = k % m_n_directions;
  int k1 = (k + 1) % m_n_directions;

  return column[k0] * (1.0 - frac) + column[k1] * frac;
}

// The per-timestep shadow test, Lambertian cosine projection (max(0, normal . sun)), and
// inverse-square distance scaling are adapted from solshade's compute_flux_timeseries
// (solshade/irradiance.py): https://github.com/amanchokshi/solshade (MIT License,
// (c) 2025 Aman Chokshi; Chokshi et al., JOSS, doi:10.21105/joss.09944). Here they are
// integrated over the diurnal cycle to a daily energy, and the solar position and
// distance factor come from PISM's own DEBMSimplePointwise instead of an ephemeris.
void TerrainInsolation::daily_insolation(double declination, double distance_factor,
                                         const array::Scalar &latitude,
                                         array::Scalar &result) const {
  const double seconds_per_day = 86400.0;

  // number of sub-daily samples used to integrate the diurnal cycle
  int M = static_cast<int>(std::lround(seconds_per_day / m_ephemeris_dt));
  if (M < 1) {
    M = 1;
  }
  const double dt = seconds_per_day / M;

  array::AccessScope scope{ &latitude, &result, m_normal_e.get(), m_normal_n.get(),
                            m_normal_u.get(), m_horizon.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    const double lat = latitude(i, j) * (M_PI / 180.0); // degrees north -> radians
    const double nE = (*m_normal_e)(i, j);
    const double nN = (*m_normal_n)(i, j);
    const double nU = (*m_normal_u)(i, j);
    const double *column = m_horizon->get_column(i, j);

    double energy = 0.0;
    for (int m = 0; m < M; ++m) {
      // hour angle sweeps the full day, midpoint rule over [-pi, pi)
      double H = -M_PI + 2.0 * M_PI * (m + 0.5) / M;

      double altitude = 0.0, azimuth = 0.0;
      terrain::sun_position(lat, declination, H, altitude, azimuth);

      if (altitude <= 0.0) {
        continue; // Sun below the astronomical horizon
      }
      if (altitude <= horizon_at(column, azimuth)) {
        continue; // shadowed by surrounding terrain
      }

      double cos_alt = std::cos(altitude);
      double sE = cos_alt * std::sin(azimuth);
      double sN = cos_alt * std::cos(azimuth);
      double sU = std::sin(altitude);

      double mu = nE * sE + nN * sN + nU * sU;
      if (mu <= 0.0) {
        continue; // surface faces away from the Sun
      }

      energy += m_solar_constant * distance_factor * mu * dt;
    }

    result(i, j) = energy;
  }
}

} // end of namespace surface
} // end of namespace pism
