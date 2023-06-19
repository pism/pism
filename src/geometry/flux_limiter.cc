/* Copyright (C) 2022, 2023 PISM Authors
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
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
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
 */
void make_nonnegative_preserving(double dt,
                                 const array::Scalar1 &x,
                                 const array::Staggered1 &flux,
                                 array::Staggered &result) {

  using details::pp;
  using details::np;
  using details::flux_out;

  auto grid = result.grid();

  double
    eps = std::numeric_limits<double>::epsilon(),
    dx = grid->dx(),
    dy = grid->dy();

  array::AccessScope list{&flux, &x, &result};

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

    const double
      div_Q   = div(Q),
      div_Q_e = div(Q_e),
      div_Q_n = div(Q_n);

    if ((div_Q   <= 0.0 or x(i, j)     - dt * div_Q   >= eps) and
        (div_Q_e <= 0.0 or x(i + 1, j) - dt * div_Q_e >= eps) and
        (div_Q_n <= 0.0 or x(i, j + 1) - dt * div_Q_n >= eps)) {
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

    // compute total amounts moved *out* of the current cell and its north and east
    // neighbors over the course of the time step dt
    //
    // see equation (A4) in [Smolarkiewicz1989]
    //
    // note that we can compute all these using the width=1 stencil because of the way
    // PISM's staggered grid is set up
    double F_out   = flux_out(Q, dx, dy, dt);
    double F_out_n = flux_out(Q_n, dx, dy, dt);
    double F_out_e = flux_out(Q_e, dx, dy, dt);

    // amounts moved through the eastern and northern cell faces
    double F_e = Q.e * dt / dx;
    double F_n = Q.n * dt / dy;

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
    double X_ij = pp(x(i, j) - eps);
    double X_e  = pp(x(i + 1, j) - eps);
    double X_n  = pp(x(i, j + 1) - eps);

    // limit total amounts (see equation (10) in [Smolarkiewicz1989])
    double F_e_limited = std::max(std::min(F_e, (pp(F_e) / F_out) * X_ij),
                                  (-np(F_e) / F_out_e) * X_e);

    assert(x(i, j) - F_e_limited >= 0);
    assert(x(i + 1, j) + F_e_limited >= 0);

    double F_n_limited = std::max(std::min(F_n, (pp(F_n) / F_out) * X_ij),
                                  (-np(F_n) / F_out_n) * X_n);

    assert(x(i, j) - F_n_limited >= 0);
    assert(x(i, j + 1) + F_n_limited >= 0);

    // convert back to fluxes:
    result(i, j, 0) = F_e_limited * dx / dt;
    result(i, j, 1) = F_n_limited * dy / dt;
  }

  limiter_count = GlobalSum(grid->com, limiter_count);
  if (limiter_count > 0) {
    auto log = grid->ctx()->log();
    log->message(2, "limited ice flux at %d locations\n", limiter_count);
  }
}

} // end of namespace pism
