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

#include "MPDATA2.hh"

#include "pism/util/array/CellType.hh"

/*! References:
 *
 * [Smolarkiewicz1983] P. K. Smolarkiewicz, “A Simple Positive Definite Advection Scheme
 * with Small Implicit Diffusion,” Monthly Weather Review, vol. 111, Art. no. 3, Mar.
 * 1983.
 *
 * [Smolarkiewicz1990] P. K. Smolarkiewicz and W. W. Grabowski, “The multidimensional
 * positive definite advection transport algorithm: nonoscillatory option,” Journal of
 * Computational Physics, vol. 86, Art. no. 2, Feb. 1990.
 */

namespace pism {

namespace fct {
// positive part
static double pp(double x) {
  return std::max(x, 0.0);
}

// negative part
static double np(double x) {
  return std::min(x, 0.0);
}

static double maximum(const stencils::Star<double> &psi) {
  using std::max;
  return max(max(max(max(psi.c, psi.n), psi.e), psi.s), psi.w);
}

static double minimum(const stencils::Star<double> &psi) {
  using std::min;
  return min(min(min(min(psi.c, psi.n), psi.e), psi.s), psi.w);
}

static double flux_in(const stencils::Star<double> &u,
                      const stencils::Star<double> &psi,
                      double dx, double dy, double dt) {
  return dt * ((pp(u.w) * psi.w - np(u.e) * psi.e) / dx +
               (pp(u.s) * psi.s - np(u.n) * psi.n) / dy);
}

static double flux_out(const stencils::Star<double> &u,
                       const stencils::Star<double> &psi,
                       double dx, double dy, double dt) {
  return dt * ((pp(u.e) * psi.c - np(u.w) * psi.c) / dx +
               (pp(u.n) * psi.c - np(u.s) * psi.c) / dy);
}

static double beta_up(const stencils::Star<double> &u,
                      const stencils::Star<double> &psi,
                      const stencils::Star<double> &psi_old,
                      double dx, double dy, double dt) {
  const double eps = 1e-15;
  double psi_max = std::max(maximum(psi), maximum(psi_old));

  return (psi_max - psi.c) / (flux_in(u, psi, dx, dy, dt) + eps);
}

static double beta_down(const stencils::Star<double> &u,
                        const stencils::Star<double> &psi,
                        const stencils::Star<double> &psi_old,
                        double dx, double dy, double dt) {
  const double eps = 1e-15;
  double psi_min = std::min(minimum(psi), minimum(psi_old));

  return (psi.c - psi_min) / (flux_out(u, psi, dx, dy, dt) + eps);
}

static void limit(double dt,
                  const array::Scalar2 &x_old,
                  const array::Scalar2 &x,
                  const array::Staggered1 &velocity,
                  array::Staggered &result) {

  auto grid = x_old.grid();

  double
    dx  = grid->dx(),
    dy  = grid->dy();

  array::AccessScope list{&x_old, &x, &velocity, &result};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double beta_u{0.0};
    double beta_d{0.0};
    {
      stencils::Star<double>
        X     = x.star(i, j),
        X_old = x_old.star(i, j),
        u     = velocity.star(i, j);

      beta_u = beta_up(u, X, X_old, dx, dy, dt);
      beta_d = beta_down(u, X, X_old, dx, dy, dt);
    }

    // east
    {
      stencils::Star<double>
        X     = x.star(i + 1, j),
        X_old = x_old.star(i + 1, j),
        u     = velocity.star(i + 1, j);

      double C{1.0};
      // note: the "west" velocity of the east neighbor is the "east" velocity of the
      // current cell
      if (u.w > 0.0) {
        C = std::min(1.0, std::min(beta_d, beta_up(u, X, X_old, dx, dy, dt)));
      } else {
        C = std::min(1.0, std::min(beta_u, beta_down(u, X, X_old, dx, dy, dt)));
      }
      result(i, j, 0) = u.w * C;
    }

    // north
    {
      stencils::Star<double>
        X     = x.star(i, j + 1),
        X_old = x_old.star(i, j + 1),
        u     = velocity.star(i, j + 1);

      double C{1.0};
      // note: the "south" velocity of the north neighbor is the "north" velocity of the
      // current cell
      if (u.s > 0.0) {
        C = std::min(1.0, std::min(beta_d, beta_up(u, X, X_old, dx, dy, dt)));
      } else {
        C = std::min(1.0, std::min(beta_u, beta_down(u, X, X_old, dx, dy, dt)));
      }
      result(i, j, 1) = u.s * C;
    }
  } // end of the loop over grid points
}

} // end of namespace fct

const array::Scalar& MPDATA2::x() const {
  return m_x;
}

MPDATA2::MPDATA2(std::shared_ptr<const Grid> grid, int N)
  : m_v(grid, "velocity_staggered"),
    m_v_old(grid, "tmp"),
    m_v_ghosted(grid, "velocity_ghosted"),
    m_x_previous(grid, "previous_state"),
    m_x_input(grid, "input_field"),
    m_x(grid, "new_state"),
    m_N(N) {
  // empty
}

/*!
 * Compute staggered grid (cell interface) velocities given regular grid velocities and
 * the domain extent (cell type).
 */
static void compute_interface_velocity(const array::CellType &cell_type,
                                       const array::Vector &velocity,
                                       array::Staggered &result) {

  auto grid = cell_type.grid();

  array::AccessScope scope{&cell_type, &velocity, &result};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    Vector2d
      V   = velocity(i, j),
      V_e = velocity(i + 1, j),
      V_n = velocity(i, j + 1);

    int
      W   = static_cast<int>(cell_type.icy(i, j)),
      W_n = static_cast<int>(cell_type.icy(i, j + 1)),
      W_e = static_cast<int>(cell_type.icy(i + 1, j));

    result(i, j, 0) = (W * V.u + W_e * V_e.u) / std::max(W + W_e, 1);
    result(i, j, 1) = (W * V.v + W_n * V_n.v) / std::max(W + W_n, 1);
  }
}

/*!
 * Compute the "anti-diffusive" flux correction in the form of velocities at cell
 * interfaces.
 */
static void compute_corrective_velocity(double dt,
                                        const array::Staggered &velocity,
                                        const array::Scalar2 &x,
                                        array::Staggered &result) {
  using std::fabs;

  auto grid = x.grid();

  const double
    eps = 1e-15,
    dx  = grid->dx(),
    dy  = grid->dy();

  array::AccessScope scope{&velocity, &x, &result};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto X = x.box(i, j);

    // eastern interface
    {
      double
        u     = velocity(i, j, 0),
        v_bar = 0.25 * (velocity(i + 1, j, 1) + velocity(i, j, 1) +
                        velocity(i, j - 1, 1) + velocity(i + 1, j - 1, 1));
      // double X_e2 = x(i + 2, j);

      // equation (13) in Smolarkiewicz1983 for the corrective velocity in the X direction
      double U = ((fabs(u) * dx - dt * u * u) * (X.e - X.c) / ((X.e + X.c + eps) * dx)
                  - 0.5 * dt * u * v_bar
                  * (X.ne - X.se + X.n - X.s) / ((X.ne + X.se + X.n + X.s + eps) * dy));

      result(i, j, 0) = U;
    }

    // northern interface
    {
      double
        v = velocity(i, j, 1),
        u_bar = 0.25 * (velocity(i, j + 1, 0) + velocity(i, j, 0) +
                        velocity(i - 1, j, 0) + velocity(i - 1, j + 1, 0));
      // double X_n2 = x(i, j + 2);

      // equation (13) in Smolarkiewicz1983 for the corrective velocity in the Y direction
      double V = ((fabs(v) * dy - dt * v * v) * (X.n - X.c) / ((X.n + X.c + eps) * dy)
                  - 0.5 * dt * v * u_bar *
                  (X.ne - X.nw + X.e - X.w) / ((X.ne + X.nw + X.e + X.w + eps) * dx));

      result(i, j, 1) = V;
    }
  }
}

/*!
 * Upwinded flux
 */
static double upwind(double x, double x_n, double u) {
  return u * (u >= 0.0 ? x : x_n);
}

/*!
 * Perform an explicit step using first order upwinding.
 *
 * @param[in] dt time step length
 * @param[in] velocity cell interface velocities
 * @param[in] x_old current state
 * @param[out] x new state
 */
static void step(double dt,
                 const array::Staggered1 &velocity,
                 const array::Scalar1 &x_old,
                 array::Scalar &x) {

  auto grid = x.grid();
  double
    dx = grid->dx(),
    dy = grid->dy();

  array::AccessScope scope{&velocity, &x_old, &x};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto u = velocity.star(i, j);
    auto f = x_old.star(i, j);

    double
      Q_e = upwind(f.c, f.e, u.e),
      Q_w = upwind(f.w, f.c, u.w),
      Q_n = upwind(f.c, f.n, u.n),
      Q_s = upwind(f.s, f.c, u.s);

    x(i, j) = x_old(i, j) - dt * ((Q_e - Q_w) / dx + (Q_n - Q_s) / dy);
  }
}

void MPDATA2::update(double dt,
                     const array::CellType &cell_type,
                     const array::Scalar &x,
                     const array::Vector &velocity,
                     bool nonoscillatory) {

  // make a ghosted copy (needed to limit fluxes)
  m_x_input.copy_from(x);

  for (int k = 0; k < m_N; ++k) {
    if (k == 0) {
      m_x_previous.copy_from(x);
      m_v_ghosted.copy_from(velocity);
      compute_interface_velocity(cell_type, m_v_ghosted, m_v);
    } else {
      m_x_previous.copy_from(m_x);
      m_v_old.copy_from(m_v);
      compute_corrective_velocity(dt, m_v_old, m_x_previous, m_v);
      if (nonoscillatory) {
        m_v_old.copy_from(m_v);
        fct::limit(dt, m_x_previous, m_x_input, m_v_old, m_v);
      }
    }

    m_v.update_ghosts();
    step(dt, m_v, m_x_previous, m_x);
  }
}

} // end of namespace pism
