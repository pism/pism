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

#include <functional>
#include <memory>
#include <vector>

#include "pism/util/petscwrappers/VecScatter.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/OrbitalParameters.hh"

namespace pism {

class Grid;

namespace array {
class Scalar;
class Scalar1;
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
  TerrainInsolation(std::shared_ptr<const Grid> grid,
                    std::function<double(double)> atmosphere_transmissivity);

  void update_horizon_map(const array::Scalar1 &surface_elevation);

  //! Daily-mean terrain-shaded surface insolation rate (W m-2) at time `time`: the
  //! diurnal insolation integral divided by the length of the day.
  void update_daily_insolation(double time, const array::Scalar1 &surface_elevation);

  const array::Scalar& insolation() const;

  void insolation_energy_series(int i, int j,
                                const std::vector<OrbitalParameters> &orbital,
                                double dt_sub, double latitude,
                                std::vector<double> &result) const;

  //! Terrain horizon map (azimuth, y, x), elevation angle in radians.
  const array::Array3D &horizon() const;

  //! Sky-view factor: the fraction (between 0 and 1) of the diffuse sky hemisphere
  //! visible from each cell, accounting for the terrain horizon and the surface
  //! slope/aspect.
  //!
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
  double m_insolation_dt;
  double m_solar_constant;
  double m_diffuse_fraction;

  petsc::VecScatter m_scatter;
  // the global DEM, replicated on every rank
  petsc::Vec m_dem_local;
  
  std::shared_ptr<array::Array3D> m_horizon;      // (azimuth, y, x), radians
  std::shared_ptr<array::Scalar> m_sky_view;      // sky-view factor, in [0, 1]

  //! daily surface insolation energy (J m-2) computed for the current update
  array::Scalar m_insolation;

  //! azimuth of the Y direction on the grid:
  array::Scalar m_y_azimuth;
  
  static double interpolate(const double *column, int n, double azimuth);

  OrbitalParameterCalculator m_orbital_parameters;

  std::function<double(double)> m_transmissivity;
};

} // end of namespace surface
} // end of namespace pism

#endif /* PISM_TERRAIN_INSOLATION_H */
