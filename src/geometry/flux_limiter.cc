/* Copyright (C) 2022, 2023, 2025 PISM Authors
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

#include "pism/geometry/flux_limiter.hh"

#include <cassert>
#include <algorithm>
#include <limits>

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Staggered.hh"
#include "pism/util/pism_utilities.hh" // GlobalSum()

namespace pism {

/*! References:
 *
 * [Smolarkiewicz1989] P. K. Smolarkiewicz, “Comment on “A Positive Definite Advection
 * Scheme Obtained by Nonlinear Renormalization of the Advective Fluxes”,” Monthly Weather
 * Review, vol. 117, Art. no. 11, 1989.
 *
 * [Zalesak1979] S. T. Zalesak, “Fully multidimensional flux-corrected transport
 * algorithms for fluids,” Journal of Computational Physics, vol. 31, Art. no. 3, Jun.
 * 1979.
 */

namespace details {

// positive part
static inline double pp(double x) {
  return std::max(x, 0.0);
}

// negative part
//
// Note the negative sign!
static inline double np(double x) {
  return -std::min(x, 0.0);
}

// Total flux out of a cell over a time step dt
static inline double flux_out(const stencils::Star<double> &u, double dx, double dy, double dt) {
  return dt * ((pp(u.e) + np(u.w)) / dx + (pp(u.n) + np(u.s)) / dy);
}

} // end of namespace details

/*! Limit fluxes to preserve non-negativity of a transported quantity.
 *
 * This method is described in [Smolarkiewicz1989].
 *
 * It is based on the [Zalesak1979] flux corrected transport limiter, but for the
 * "regular" flux instead of the "anti-diffusive" flux and with a different limiting
 * criterion (non-negativity instead of monotonicity).
 *
 * @param[in] Q_c fluxes through sides of the current ("c") cell
 * @param[in] Q_e fluxes through sides of the eastern ("e") neighbor of the current cell
 * @param[in] Q_n fluxes through sides of the northern ("n") neighbor of the current cell
 * @param[in] x_c value at the current cell
 * @param[in] x_e value at the eastern neighbor
 * @param[in] x_n value at the northern neighbor
 * @param[in] dx grid spacing in the X direction
 * @param[in] dy grid spacing in the Y direction
 * @param[in] dt time step length
 * @param[in] eps lower bound of the transported quantity (a small positive constant)
 *
 * Returns {Q_e, Q_n} - limited fluxed through eastern and northern sides of the current
 * grid cell.
 */
std::array<double, 2> flux_limiter(const stencils::Star<double> &Q_c,
                                   const stencils::Star<double> &Q_e,
                                   const stencils::Star<double> &Q_n, double x_c, double x_e,
                                   double x_n, double dx, double dy, double dt, double eps) {

  using details::flux_out;
  using details::np;
  using details::pp;

  // compute total amounts moved *out* of the current cell and its north and east
  // neighbors over the course of the time step dt
  //
  // see equation (A4) in [Smolarkiewicz1989]
  //
  // note that we can compute all these using the width=1 stencil because of the way
  // PISM's staggered grid is set up
  double F_out   = flux_out(Q_c, dx, dy, dt);
  double F_out_n = flux_out(Q_n, dx, dy, dt);
  double F_out_e = flux_out(Q_e, dx, dy, dt);

  // amounts moved through the eastern and northern cell faces
  double F_e = Q_c.e * dt / dx;
  double F_n = Q_c.n * dt / dy;

  // Maximum amounts the current cell and its neighbors can lose while maintaining
  // non-negativity
  //
  // Note: we limit total amounts so that
  //
  // - if a cell value X is below eps, the flux is zero
  //
  // - otherwise the total flux out of a cell can remove at most (X - eps) over the
  //   course of a time step
  //
  // This is needed to avoid small negative values resulting from rounding errors.
  double X_c = pp(x_c - eps);
  double X_e = pp(x_e - eps);
  double X_n = pp(x_n - eps);

  // limit total amounts (see equation (10) in [Smolarkiewicz1989])
  double F_e_limited = std::max(std::min(F_e, (pp(F_e) / F_out) * X_c), (-np(F_e) / F_out_e) * X_e);

  assert(x_c - F_e_limited >= 0);
  assert(x_e + F_e_limited >= 0);

  double F_n_limited = std::max(std::min(F_n, (pp(F_n) / F_out) * X_c), (-np(F_n) / F_out_n) * X_n);

  assert(x_c - F_n_limited >= 0);
  assert(x_n + F_n_limited >= 0);

  // convert back to fluxes and return:
  return { F_e_limited * dx / dt, F_n_limited * dy / dt };
}

/*! Limit fluxes to preserve non-negativity of a transported quantity.
 *
 * See flux_limiter() for details.
 */
int make_nonnegative_preserving(double dt, const array::Scalar1 &x, const array::Staggered1 &flux,
                                array::Staggered &result) {

  auto grid = result.grid();

  double eps = std::numeric_limits<double>::epsilon(), dx = grid->dx(), dy = grid->dy();

  array::AccessScope list{ &flux, &x, &result };

  // flux divergence
  auto div = [dx, dy](const stencils::Star<double> &Q) {
    return (Q.e - Q.w) / dx + (Q.n - Q.s) / dy;
  };

  int limiter_count = 0;

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto Q   = flux.star(i, j);
    auto Q_n = flux.star(i, j + 1);
    auto Q_e = flux.star(i + 1, j);

    double x_c = x(i, j);
    double x_e = x.E(i, j);
    double x_n = x.N(i, j);

    const double div_Q = div(Q), div_Q_e = div(Q_e), div_Q_n = div(Q_n);

    if ((div_Q <= 0.0 or x_c - dt * div_Q >= eps) and
        (div_Q_e <= 0.0 or x_e - dt * div_Q_e >= eps) and
        (div_Q_n <= 0.0 or x_n - dt * div_Q_n >= eps)) {
      // No need to limit fluxes: total fluxes out of cells (i, j), (i + 1, j), (i, j + 1)
      // may be able to create a negative thickness, but fluxes *into* these cells make up for it
      //
      // Without this check the limiter is a little over-eager and may affect results in
      // areas where mass conservation is not an issue.
      result(i, j, 0) = Q.e;
      result(i, j, 1) = Q.n;
      continue;
    }

    limiter_count += 1;

    auto Q_l = flux_limiter(Q, Q_e, Q_n, x_c, x_e, x_n, dx, dy, dt, eps);

    result(i, j, 0) = Q_l[0];
    result(i, j, 1) = Q_l[1];
  }

  return GlobalSum(grid->com, limiter_count);
}

} // end of namespace pism
