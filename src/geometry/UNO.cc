/* Copyright (C) 2020, 2022, 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <cmath>

#include "UNO.hh"
#include "flux_limiter.hh"

#include "pism/util/array/CellType.hh"

/*! References:
 *
 * [Li2008] J.-G. Li, “Upstream Nonoscillatory Advection Schemes,” Monthly Weather Review,
 * vol. 136, Art. no. 12, Dec. 2008.
 */

namespace pism {

// 1D Upwind Non-Oscillatory (UNO) advection schemes and related approximations
namespace uno {

// the difference quotient (equation 2 in [Li2008])
static inline double difference_quotient(const double *x, const double *f, size_t A, size_t B) {
  return (f[A] - f[B]) / (x[A] - x[B]);
}

// index of the "upstream" cell
static inline size_t upstream(double v, size_t j) {
  return v > 0.0 ? j - 1 : j + 2;
}

// index of the "central" cell
static inline size_t central(double v, size_t j) {
  return v > 0.0 ? j : j + 1;
}

// index of the "downstream" cell
static inline size_t downstream(double v, size_t j) {
  return v > 0.0 ? j + 1 : j;
}

// the sign
static inline double sign(double x) {
  return x >= 0.0 ? 1.0 : -1.0;
}

// first order upwinding
static inline double psi_upwind1(const double *x,
                                 const double *y,
                                 size_t j,
                                 double velocity,
                                 double dx,
                                 double dt) {
  (void) x;
  (void) dx;
  (void) dt;

  return y[central(velocity, j)];
}

// Lax-Wendroff
static inline double psi_lax_wendroff(const double *x,
                                      const double *y,
                                      size_t j,
                                      double velocity,
                                      double dx,
                                      double dt) {
  size_t
    c = central(velocity, j),
    d = downstream(velocity, j);

  double gc = difference_quotient(x, y, d, c);

  // equation 3 in [Li2008]
  return y[c] + 0.5 * sign(velocity) * (dx - std::abs(velocity) * dt) * gc;
}

// Fromm
static inline double psi_fromm(const double *x,
                               const double *y,
                               size_t j,
                               double velocity,
                               double dx,
                               double dt) {
  size_t
    u = upstream(velocity, j),
    c = central(velocity, j),
    d = downstream(velocity, j);

  double gc = difference_quotient(x, y, d, u);

  // equation 3 in [Li2008]
  return y[c] + 0.5 * sign(velocity) * (dx - std::abs(velocity) * dt) * gc;
}

// UNO2
static inline double psi_uno2(const double *x,
                              const double *y,
                              size_t j,
                              double velocity,
                              double dx,
                              double dt) {
  size_t
    u = upstream(velocity, j),
    c = central(velocity, j),
    d = downstream(velocity, j);

  double
    gdc = difference_quotient(x, y, d, c),
    gcu = difference_quotient(x, y, c, u);

  // equation 5 in [Li2008]
  double gc = sign(gdc) * std::min(std::abs(gdc), std::abs(gcu));

  // equation 3 in [Li2008]
  return y[c] + 0.5 * sign(velocity) * (dx - std::abs(velocity) * dt) * gc;
}

// UNO3
static inline double psi_uno3(const double *x,
                              const double *y,
                              size_t j,
                              double velocity,
                              double dx,
                              double dt) {
  size_t
    u = upstream(velocity, j),
    c = central(velocity, j),
    d = downstream(velocity, j);

  double
    gdc = difference_quotient(x, y, d, c),
    gcu = difference_quotient(x, y, c, u),
    gdu = difference_quotient(x, y, d, u);

  double abs_u = std::abs(velocity);
  double sgn_u = sign(velocity);

  // equation 10 in [Li2008]
  double gc = 0.0;
  // NOLINTBEGIN(readability-magic-numbers)
  if (std::abs(gdc - gcu) < 1.2 * std::abs(gdu)) {
    gc = gdc - ((dx + abs_u * dt) / (1.5 * sgn_u)) * ((gdc - gcu) / (x[d] - x[u]));
  } else if (gdc * gcu > 0.0) {
    gc =  2 * sign(gdc) * std::min(std::abs(gdc), std::abs(gcu));
  } else {
    gc = sign(gdc) * std::min(std::abs(gdc), std::abs(gcu)); // UNO2
  }
  // NOLINTEND(readability-magic-numbers)

  return y[c] + 0.5 * sgn_u * (dx - abs_u * dt) * gc;
}

} // end of namespace uno

const array::Scalar& UNO::x() const {
  return m_x;
}

UNO::UNO(std::shared_ptr<const IceGrid> grid, UNOType type)
  : m_q(grid, "interface_fluxes"),
    m_q_limited(grid, "limited_interface_fluxes"),
    m_v_ghosted(grid, "velocity"),
    m_x_ghosted(grid, "old_state"),
    m_x(grid, "new_state")
{
  switch (type) {
  case PISM_UNO_UPWIND1:
    m_approx = uno::psi_upwind1;
    break;
  case PISM_UNO_LAX_WENDROFF:
    m_approx = uno::psi_lax_wendroff;
    break;
  case PISM_UNO_FROMM:
    m_approx = uno::psi_fromm;
    break;
  case PISM_UNO_2:
    m_approx = uno::psi_uno2;
    break;
  default:
  case PISM_UNO_3:
    m_approx = uno::psi_uno3;
    break;
  }
}

/*!
 * Compute staggered grid (cell interface) fluxes given the domain extent (cell type),
 * regular grid velocities, and the current distribution of the advected quantity.
 */
void UNO::compute_interface_fluxes(const array::CellType1 &cell_type,
                                   const array::Vector1 &velocity,
                                   const array::Scalar2 &x_old,
                                   double dt,
                                   array::Staggered &result) const {

  auto grid = cell_type.grid();

  double
    dx = grid->dx(),
    dy = grid->dy();

  array::AccessScope scope{&cell_type, &velocity, &x_old, &result};

  // temporary storage for values needed by stencil computations
  double coords[4];
  double values[4];

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    // velocities through east and north cell interfaces:
    double vx, vy;
    {
      Vector2d
        V   = velocity(i, j),
        V_e = velocity(i + 1, j),
        V_n = velocity(i, j + 1);

      double
        W   = static_cast<double>(cell_type.icy(i, j)),
        W_n = static_cast<double>(cell_type.icy(i, j + 1)),
        W_e = static_cast<double>(cell_type.icy(i + 1, j));

      vx = (W * V.u + W_e * V_e.u) / std::max(W + W_e, 1.0);
      vy = (W * V.v + W_n * V_n.v) / std::max(W + W_n, 1.0);
    }

    // east
    {
      // gather inputs to the 1D UNO code: we need four numbers centered around i+1/2,
      // i.e. [i-1, i, i+1, i+2].
      double x0 = i > 0 ? grid->x(i - 1) : grid->x(0) - dx;
      for (int k = 0; k < 4; ++k) {
        coords[k] = x0 + static_cast<double>(k) * dx;
        values[k] = x_old((i - 1) + k, j);
      }

      // use the 1D "mid-flux" approximation to get the flux
      result(i, j, 0) = vx * m_approx(coords, values, 1, vx, dx, dt);
    }

    // north
    {
      // gather inputs to the 1D UNO code: we need four numbers centered around j+1/2,
      // i.e. [j-1, j, j+1, j+2].
      double y0 = j > 0 ? grid->y(j - 1) : grid->y(0) - dy;
      for (int k = 0; k < 4; ++k) {
        coords[k] = y0 + k * dy;
        values[k] = x_old(i, (j - 1) + k);
      }

      // use the 1D "mid-flux" approximation to get the flux
      result(i, j, 1) = vy * m_approx(coords, values, 1, vy, dy, dt);
    }
  } // end of the loop over grid points
}

/*!
 * Perform an explicit step using the flux form of the advection equation.
 *
 * @param[in] dt time step length
 * @param[in] flux cell interface fluxes
 * @param[in] x_old current state
 * @param[out] result new state
 */
static void step(double dt,
                 const array::Staggered1 &flux,
                 const array::Scalar &x_old,
                 array::Scalar &result) {

  auto grid = result.grid();

  double
    dx  = grid->dx(),
    dy  = grid->dy();

  array::AccessScope scope{&flux, &x_old, &result};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const auto Q = flux.star(i, j);

    result(i, j) = x_old(i, j) - dt * ((Q.e - Q.w) / dx + (Q.n - Q.s) / dy);
  }
}

void UNO::update(double dt,
                 const array::CellType1 &cell_type,
                 const array::Scalar &x,
                 const array::Vector &velocity,
                 bool nonnegative) {

  // make ghosted copies:
  m_v_ghosted.copy_from(velocity);
  m_x_ghosted.copy_from(x);

  compute_interface_fluxes(cell_type, m_v_ghosted, m_x_ghosted, dt, m_q);

  m_q.update_ghosts();

  // limit fluxes to preserve non-negativity
  if (nonnegative) {
    make_nonnegative_preserving(dt, m_x_ghosted, m_q, m_q_limited);
    m_q.copy_from(m_q_limited);
  }

  step(dt, m_q, x, m_x);
}

} // end of namespace pism
