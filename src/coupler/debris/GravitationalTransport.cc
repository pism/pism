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

#include <algorithm>
#include <cmath>
#include <limits>

#include "pism/coupler/debris/GravitationalTransport.hh"

#include "pism/geometry/flux_limiter.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace debris {

GravitationalTransport::GravitationalTransport(std::shared_ptr<const Grid> grid, double mobility,
                                               double solid_density, double cfl_ratio)
  : m_grid(grid),
    m_flux(grid, "debris_gravitational_flux_staggered"),
    m_flux_limited(grid, "debris_gravitational_flux_limited"),
    m_surface(grid, "debris_surface_elevation"),
    m_cfl_ratio(cfl_ratio),
    m_substeps(0) {

  const double g = grid->ctx()->config()->get_number("constants.standard_gravity");

  m_coefficient = mobility * solid_density * g;

  m_flux.set(0.0);
  m_flux_limited.set(0.0);
}

double GravitationalTransport::coefficient() const {
  return m_coefficient;
}

double GravitationalTransport::max_dt(double K_max) const {
  if (K_max <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  const double dx = m_grid->dx(), dy = m_grid->dy();
  return m_cfl_ratio / (K_max * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
}

/*!
 * Flux across the eastern and northern faces of every cell: `-K grad(phi)`, where
 * `phi = h_s + h_d` and `K = coefficient * h_d` is taken from the cell the flux leaves
 * (upwind). This makes the flux out of a cell without debris exactly zero, so debris does
 * not "spread" from empty cells (which would also break mass conservation because the
 * non-negativity limiter assumes that an inflow it sees will actually arrive). Faces to
 * ice-free cells are closed.
 *
 * @return the largest "velocity" `coefficient * |grad(phi)|` across a face (for the
 * advective stability limit)
 */
double GravitationalTransport::compute_flux(const array::CellType1 &cell_type,
                                            const array::Scalar1 &h_d) {
  const double dx = m_grid->dx(), dy = m_grid->dy();

  array::AccessScope list{ &cell_type, &h_d, &m_surface, &m_flux };

  double speed_max = 0.0;

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    const bool icy = cell_type.icy(i, j);
    const double phi = m_surface(i, j);

    // east
    if (icy and cell_type.icy(i + 1, j)) {
      const double slope = (m_surface(i + 1, j) - phi) / dx;
      // the flux goes down the slope: from c to e if slope < 0
      const double K = m_coefficient * (slope < 0.0 ? h_d(i, j) : h_d(i + 1, j));
      m_flux(i, j, 0) = -K * slope;
      speed_max = std::max(speed_max, m_coefficient * std::fabs(slope));
    } else {
      m_flux(i, j, 0) = 0.0;
    }

    // north
    if (icy and cell_type.icy(i, j + 1)) {
      const double slope = (m_surface(i, j + 1) - phi) / dy;
      const double K = m_coefficient * (slope < 0.0 ? h_d(i, j) : h_d(i, j + 1));
      m_flux(i, j, 1) = -K * slope;
      speed_max = std::max(speed_max, m_coefficient * std::fabs(slope));
    } else {
      m_flux(i, j, 1) = 0.0;
    }
  }

  return GlobalMax(m_grid->com, speed_max);
}

void GravitationalTransport::update_surface(const array::Scalar &surface_elevation,
                                            const array::Scalar1 &h_d) {
  array::AccessScope list{ &surface_elevation, &h_d, &m_surface };
  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    m_surface(i, j) = surface_elevation(i, j) + h_d(i, j);
  }
  m_surface.update_ghosts();
}

void GravitationalTransport::step(double dt, const array::CellType1 &cell_type,
                                  const array::Scalar &surface_elevation, array::Scalar1 &h_d) {
  m_substeps = 0;

  const double K_max = m_coefficient * array::max(h_d);
  if (K_max <= 0.0 or dt <= 0.0) {
    m_flux.set(0.0);
    m_flux_limited.set(0.0);
    return;
  }

  const double dx = m_grid->dx(), dy = m_grid->dy();

  // Stability: the diffusive limit from K and, since the slope of the ice surface makes
  // this an advection with the speed coefficient * |grad(h_s)| where the debris is thin,
  // the advective (CFL) limit from the initial slope of the debris surface.
  update_surface(surface_elevation, h_d);
  const double speed_max = compute_flux(cell_type, h_d);
  double dt_max = max_dt(K_max);
  if (speed_max > 0.0) {
    dt_max = std::min(dt_max, m_cfl_ratio * std::min(dx, dy) / speed_max);
  }

  const int n = std::max(1, (int)std::ceil(dt / dt_max));
  const double dt_sub = dt / n;

  for (int s = 0; s < n; ++s) {
    if (s > 0) {
      update_surface(surface_elevation, h_d);
      compute_flux(cell_type, h_d);
    }
    m_flux.update_ghosts();

    make_nonnegative_preserving(dt_sub, h_d, m_flux, m_flux_limited);
    m_flux_limited.update_ghosts();

    {
      array::AccessScope list{ &m_flux_limited, &h_d };
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();
        auto Q = m_flux_limited.star(i, j);
        h_d(i, j) -= dt_sub * ((Q.e - Q.w) / dx + (Q.n - Q.s) / dy);
        // rounding can leave tiny negative values behind
        h_d(i, j) = std::max(h_d(i, j), 0.0);
      }
    }
    h_d.update_ghosts();

    m_substeps += 1;
  }
}

const array::Staggered1 &GravitationalTransport::flux() const {
  return m_flux_limited;
}

int GravitationalTransport::substeps() const {
  return m_substeps;
}

} // end of namespace debris
} // end of namespace pism
