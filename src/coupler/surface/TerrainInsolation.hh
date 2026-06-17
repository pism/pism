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

#ifndef PISM_TERRAIN_INSOLATION_H
#define PISM_TERRAIN_INSOLATION_H

#include <memory>
#include <vector>

namespace pism {

class Grid;

namespace array {
class Scalar;
class Array3D;
} // namespace array

namespace surface {

/*!
 * Computes the terrain horizon map and terrain-shaded daily surface insolation from a
 * digital elevation model (DEM), porting the Python "solshade" algorithms into PISM. Used
 * by the `debm_enhanced` surface model when no precomputed insolation file is given.
 *
 * The terrain horizon ray-casting reaches several kilometres from each cell, crossing MPI
 * subdomain boundaries. Because the DEM is a single, small 2D field, `init()` gathers the
 * full global DEM onto every rank (proc0 gather + broadcast); each rank then ray-marches
 * only its owned cells against the shared global DEM with no ghost communication.
 *
 * The horizon (azimuth, y, x) and the surface normals are computed once in `init()`. The
 * daily insolation field is computed on demand by `daily_insolation()` from the analytic
 * solar declination and Sun-Earth distance factor for a given day; the daily integral over
 * the diurnal cycle is independent of longitude, so only latitude, declination, distance
 * factor, and the per-cell horizon and normal are needed.
 */
class TerrainInsolation {
public:
  TerrainInsolation(std::shared_ptr<const Grid> grid);

  //! Gather the DEM and compute the horizon map and surface normals (one-time).
  void init(const array::Scalar &surface_elevation);

  //! Daily-mean terrain-shaded surface insolation rate (W m-2) for a day with the given
  //! solar `declination` (radians) and `distance_factor` (= (d_bar/d)^2): the diurnal
  //! insolation integral divided by the length of the day. `latitude` is the per-cell
  //! latitude field (degrees north); `result` is overwritten.
  void daily_insolation(double declination, double distance_factor,
                        const array::Scalar &latitude, array::Scalar &result) const;

  //! Terrain horizon map (azimuth, y, x), elevation angle in radians.
  const array::Array3D &horizon() const;

  //! Sky-view factor (y, x), in [0, 1]: the fraction of the diffuse sky hemisphere visible
  //! from each cell, accounting for the terrain horizon and the surface slope/aspect.
  //! Only valid (computed) when sky_view_enabled() is true.
  const array::Scalar &sky_view() const;

  //! Whether the sky-view factor is computed (surface.debm_enhanced.use_sky_view_factor).
  bool sky_view_enabled() const;

private:
  std::shared_ptr<const Grid> m_grid;

  // configuration
  int m_n_directions;
  double m_max_distance;
  double m_step;
  double m_ephemeris_dt;
  double m_solar_constant;
  bool m_use_sky_view;
  double m_diffuse_fraction;

  // azimuth sample directions (radians, clockwise from north), used as Array3D levels
  std::vector<double> m_azimuth;

  // the global DEM, replicated on every rank, row-major dem[j * Mx + i]
  std::vector<double> m_dem;

  std::shared_ptr<array::Array3D> m_horizon;      // (azimuth, y, x), radians
  std::shared_ptr<array::Scalar> m_normal_e;      // east component of surface normal
  std::shared_ptr<array::Scalar> m_normal_n;      // north component
  std::shared_ptr<array::Scalar> m_normal_u;      // up component
  std::shared_ptr<array::Scalar> m_sky_view;      // sky-view factor, in [0, 1]

  double horizon_at(const double *column, double azimuth) const;
};

} // end of namespace surface
} // end of namespace pism

#endif /* PISM_TERRAIN_INSOLATION_H */
