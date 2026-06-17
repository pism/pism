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

#ifndef PISM_TERRAIN_INSOLATION_KERNEL_H
#define PISM_TERRAIN_INSOLATION_KERNEL_H

namespace pism {
namespace surface {
namespace terrain {

/*!
 * Pure (PISM-free) kernels for the terrain-horizon and terrain-shaded insolation
 * computation ported from the Python "solshade" package. See
 * `examples/debm_enhanced/solstis.md` for the underlying math.
 *
 * Conventions used throughout:
 *   - The global DEM is row-major with `i` the x-index (increasing eastward) and `j`
 *     the y-index (increasing northward): `dem[j * Mx + i]`. This matches PISM's grid,
 *     where `grid->x(i)` increases with `i` and `grid->y(j)` increases with `j`.
 *   - Azimuth is measured clockwise from north (0 = north, pi/2 = east), so the
 *     horizontal unit vector in a given azimuth direction is (East, North) =
 *     (sin(azimuth), cos(azimuth)). The same convention is used for the sun azimuth,
 *     making the shadow test self-consistent.
 *   - All angles are in radians.
 */

//! Bilinear sample of the global DEM at fractional grid coordinates (fi, fj).
//! Coordinates are clamped to [0, Mx-1] x [0, My-1].
double sample_bilinear(const double *dem, int Mx, int My, double fi, double fj);

//! Terrain horizon elevation angle (radians) at cell (i0, j0) looking in direction
//! `azimuth` (radians, clockwise from north). Rays march in `step` (m) increments out to
//! `max_distance` (m); `dx`, `dy` are the grid spacings (m). The horizon angle is the
//! maximum of atan2(z(d) - z0, d) over the ray. Rays stop at the domain boundary. Returns
//! 0 if the ray immediately leaves the domain (cell on the boundary).
double ray_horizon(const double *dem, int Mx, int My, double dx, double dy,
                   int i0, int j0, double azimuth, double step, double max_distance);

//! Upward unit surface normal in East-North-Up components from the surface-elevation
//! gradients `dzdE = dz/dEast` and `dzdN = dz/dNorth` (both m/m).
void surface_normal(double dzdE, double dzdN, double &nE, double &nN, double &nU);

//! Sun altitude and azimuth (radians) from `latitude`, solar `declination`, and
//! `hour_angle` (all radians; hour angle measured from solar noon, positive in the
//! afternoon). `altitude` is in [-pi/2, pi/2]; `azimuth` is clockwise from north in
//! [0, 2*pi).
void sun_position(double latitude, double declination, double hour_angle,
                  double &altitude, double &azimuth);

//! Sky-view factor: the fraction of isotropic diffuse sky irradiance received by a tilted,
//! horizon-obstructed surface relative to an unobstructed horizontal surface, in [0, 1].
//!
//! Implements the slope-corrected formula of Dozier & Frew (1990). `horizon[k]` is the
//! terrain horizon elevation angle (radians, from horizontal) in azimuth `azimuth[k]`
//! (radians, clockwise from north), for `n_dir` equally-spaced directions. `slope` and
//! `aspect` (radians; aspect = downslope azimuth, clockwise from north) describe the
//! surface tilt. A horizontal, unobstructed surface returns 1; a flat surface with uniform
//! horizon elevation h returns cos(h)^2.
double sky_view_factor(const double *horizon, const double *azimuth, int n_dir,
                       double slope, double aspect);

} // end of namespace terrain
} // end of namespace surface
} // end of namespace pism

#endif /* PISM_TERRAIN_INSOLATION_KERNEL_H */
