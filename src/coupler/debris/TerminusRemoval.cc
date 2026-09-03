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

#include "pism/coupler/debris/TerminusRemoval.hh"

#include "pism/coupler/debris/terminus_kernels.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace debris {

TerminusRemoval::TerminusRemoval(std::shared_ptr<const Grid> grid, double gamma, double mobility,
                                 double solid_density)
  : m_grid(grid),
    m_surface(grid, "debris_surface_elevation_wide"),
    m_removal_rate(grid, "debris_removal_rate"),
    m_gamma(gamma) {

  const double g = grid->ctx()->config()->get_number("constants.standard_gravity");
  m_coefficient = mobility * solid_density * g;

  if (gamma <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "debris.transport.marginal_length_scale has to be positive (got %f)",
                                  gamma);
  }

  m_n = std::max(0, (int)std::ceil(gamma / grid->dx()) - 1);

  // the walk up-glacier is limited by the width of the ghost zone of m_surface and h_d
  const int max_n = 2;
  if (m_n > max_n) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "debris.transport.marginal_length_scale = %f m would average over %d cells"
                                  " up-glacier of the margin; at most %d are supported (%f m at this resolution)",
                                  gamma, m_n, max_n, (max_n + 1) * grid->dx());
  }

  m_removal_rate.metadata(0)
      .long_name("rate of removal of supraglacial debris into the foreland")
      .units("m s^-1")
      .output_units("m year^-1");
  m_removal_rate.set(0.0);
}

int TerminusRemoval::upstream_cells() const {
  return m_n;
}

void TerminusRemoval::step(double dt, const array::CellType2 &cell_type,
                           const array::Scalar &surface_elevation, array::Scalar2 &h_d) {
  const double dx = m_grid->dx(), dy = m_grid->dy();

  // surface of the debris layer, with wide ghosts for the walk up-glacier
  {
    array::AccessScope list{ &surface_elevation, &h_d, &m_surface };
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      m_surface(i, j) = surface_elevation(i, j) + h_d(i, j);
    }
  }
  m_surface.update_ghosts();

  {
    array::AccessScope list{ &cell_type, &h_d, &m_surface, &m_removal_rate };

    auto surface  = [&](int i, int j) -> double { return m_surface(i, j); };
    auto icy      = [&](int i, int j) -> bool { return cell_type.icy(i, j); };
    auto ice_free = [&](int i, int j) -> bool { return cell_type.ice_free(i, j); };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      m_removal_rate(i, j) = 0.0;

      if (not cell_type.icy(i, j) or h_d(i, j) <= 0.0) {
        continue;
      }

      const double slope = terminus::foreland_slope(surface, ice_free, i, j, dx, dy);
      if (slope <= 0.0) {
        continue;
      }

      // average the debris thickness and the slope along the steepest ascent
      double h_sum = h_d(i, j), s_sum = slope;
      int count = 1;

      int pi = i, pj = j;
      for (const auto &o : terminus::upstream_chain(surface, icy, i, j, m_n, dx, dy)) {
        const int ci = i + o.di, cj = j + o.dj;
        const double distance = std::sqrt(std::pow((ci - pi) * dx, 2) + std::pow((cj - pj) * dy, 2));

        h_sum += h_d(ci, cj);
        s_sum += std::fabs(m_surface(ci, cj) - m_surface(pi, pj)) / distance;
        count += 1;

        pi = ci;
        pj = cj;
      }

      const double h_eff = h_sum / count, s_eff = s_sum / count;

      // equation 23: flux (m^2 s^-1) over the marginal length scale is a thickness loss rate
      const double rate = m_coefficient * h_eff * s_eff / m_gamma;

      m_removal_rate(i, j) = std::min(rate, h_d(i, j) / dt);
    }
  }

  {
    array::AccessScope list{ &h_d, &m_removal_rate };
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      h_d(i, j) = std::max(h_d(i, j) - m_removal_rate(i, j) * dt, 0.0);
    }
  }
  h_d.update_ghosts();
}

const array::Scalar &TerminusRemoval::removal_rate() const {
  return m_removal_rate;
}

} // end of namespace debris
} // end of namespace pism
