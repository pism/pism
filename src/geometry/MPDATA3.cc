/* Copyright (C) 2026 PISM Authors
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

#include <algorithm>
#include <cmath>
#include <limits>

#include "pism/geometry/MPDATA3.hh"

#include "pism/geometry/mpdata_helpers.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"

/*! References:
 *
 * [Smolarkiewicz1983] P. K. Smolarkiewicz, "A Simple Positive Definite Advection Scheme
 * with Small Implicit Diffusion," Monthly Weather Review, vol. 111, no. 3, 1983.
 *
 * [Smolarkiewicz1990] P. K. Smolarkiewicz and W. W. Grabowski, "The multidimensional
 * positive definite advection transport algorithm: nonoscillatory option," Journal of
 * Computational Physics, vol. 86, no. 2, 1990.
 *
 * See MPDATA2.cc for the two-dimensional version this follows.
 */

namespace pism {

using mpdata::np;
using mpdata::pp;
using mpdata::upwind;

static const double eps = 1e-15;

MPDATA3::MPDATA3(std::shared_ptr<const Grid> grid, int N, bool nonoscillatory)
  : m_u_face(grid, "u_face", array::WITH_GHOSTS, grid->z(), 1),
    m_v_face(grid, "v_face", array::WITH_GHOSTS, grid->z(), 1),
    m_w_face(grid, "w_face", array::WITH_GHOSTS, grid->z(), 1),
    m_u_old(grid, "u_face_old", array::WITH_GHOSTS, grid->z(), 1),
    m_v_old(grid, "v_face_old", array::WITH_GHOSTS, grid->z(), 1),
    m_w_old(grid, "w_face_old", array::WITH_GHOSTS, grid->z(), 1),
    m_x_previous(grid, "previous_state", array::WITH_GHOSTS, grid->z(), 2),
    m_x_input(grid, "input_field", array::WITH_GHOSTS, grid->z(), 2),
    m_x(grid, "new_state", array::WITHOUT_GHOSTS, grid->z()),
    m_N(N),
    m_nonoscillatory(nonoscillatory) {

  if (N < 1) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "MPDATA3 requires at least one pass (got %d)", N);
  }

  // control volume interfaces (see debris::column::interfaces())
  const auto &z = grid->z();
  const size_t Mz = z.size();

  m_zi.resize(Mz + 1);
  m_zi[0] = 0.0;
  for (size_t k = 1; k < Mz; ++k) {
    m_zi[k] = 0.5 * (z[k - 1] + z[k]);
  }
  m_zi[Mz] = std::numeric_limits<double>::infinity();

  m_dz_face.resize(Mz > 0 ? Mz - 1 : 0);
  for (size_t k = 0; k + 1 < Mz; ++k) {
    m_dz_face[k] = z[k + 1] - z[k];
  }
}

const array::Array3D &MPDATA3::x() const {
  return m_x;
}

double MPDATA3::dz(int k, double H) const {
  return std::max(0.0, std::min(m_zi[k + 1], H) - m_zi[k]);
}

/*!
 * Interface velocities of the first (upwind) pass.
 *
 * A horizontal face is open only if both columns are ice-covered and contain ice at this
 * level; a vertical face is open only if both volumes contain ice. Closed faces have zero
 * velocity and stay closed in all passes because every term of the anti-diffusive
 * velocity is proportional to the velocity of the previous pass.
 */
void MPDATA3::compute_interface_velocity(const array::Scalar1 &ice_thickness,
                                         const array::CellType1 &cell_type,
                                         const array::Array3D &u, const array::Array3D &v,
                                         const array::Array3D &w) {
  auto grid = m_x.grid();
  const int Mz = (int)grid->Mz();

  array::AccessScope list{ &ice_thickness, &cell_type, &u, &v, &w,
                           &m_u_face, &m_v_face, &m_w_face };

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    const double H_c = ice_thickness(i, j), H_e = ice_thickness(i + 1, j),
                 H_n = ice_thickness(i, j + 1);
    const bool icy_c = cell_type.icy(i, j), icy_e = cell_type.icy(i + 1, j),
               icy_n = cell_type.icy(i, j + 1);

    const double *u_c = u.get_column(i, j), *u_e = u.get_column(i + 1, j),
                 *v_c = v.get_column(i, j), *v_n = v.get_column(i, j + 1),
                 *w_c = w.get_column(i, j);

    double *uf = m_u_face.get_column(i, j), *vf = m_v_face.get_column(i, j),
           *wf = m_w_face.get_column(i, j);

    for (int k = 0; k < Mz; ++k) {
      const bool open_c = icy_c and dz(k, H_c) > 0.0;

      uf[k] = (open_c and icy_e and dz(k, H_e) > 0.0) ? 0.5 * (u_c[k] + u_e[k]) : 0.0;
      vf[k] = (open_c and icy_n and dz(k, H_n) > 0.0) ? 0.5 * (v_c[k] + v_n[k]) : 0.0;

      if (k + 1 < Mz and open_c and dz(k + 1, H_c) > 0.0) {
        wf[k] = 0.5 * (w_c[k] + w_c[k + 1]);
      } else {
        wf[k] = 0.0;
      }
    }
  }
}

/*!
 * Anti-diffusive velocities (equation 13 in [Smolarkiewicz1983]) computed from the
 * velocities of the previous pass (`m_*_old`) and its result (`m_x_previous`).
 *
 * Horizontal faces use mass per unit area, vertical faces use the concentration; see the
 * class documentation.
 */
void MPDATA3::compute_corrective_velocity(double dt, const array::Scalar1 &ice_thickness) {
  auto grid = m_x.grid();
  const int Mz = (int)grid->Mz();
  const double dx = grid->dx(), dy = grid->dy();
  const auto &z = grid->z();

  array::AccessScope list{ &ice_thickness, &m_x_previous,
                           &m_u_old, &m_v_old, &m_w_old,
                           &m_u_face, &m_v_face, &m_w_face };

  // concentration in volume k of the column at (i, j) (zero where there is no ice)
  auto C = [&](int i, int j, int k) -> double {
    double t = dz(k, ice_thickness(i, j));
    return t > 0.0 ? m_x_previous.get_column(i, j)[k] / t : 0.0;
  };

  auto X = [&](int i, int j, int k) -> double {
    return m_x_previous.get_column(i, j)[k];
  };

  auto w_old = [&](int i, int j, int k) -> double {
    return k >= 0 ? m_w_old.get_column(i, j)[k] : 0.0;
  };

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    double *uf = m_u_face.get_column(i, j), *vf = m_v_face.get_column(i, j),
           *wf = m_w_face.get_column(i, j);

    for (int k = 0; k < Mz; ++k) {
      const int kp = std::min(k + 1, Mz - 1), km = std::max(k - 1, 0);
      const double dz_c = 0.5 * (z[kp] - z[km]); // half the distance between the volumes above and below

      // eastern face
      {
        const double u = m_u_old.get_column(i, j)[k];
        double U = 0.0;
        if (u != 0.0) {
          U = mpdata::corrective_velocity(dt, u, X(i, j, k), X(i + 1, j, k), dx, eps);

          const double v_bar = 0.25 * (m_v_old.get_column(i + 1, j)[k] + m_v_old.get_column(i, j)[k] +
                                       m_v_old.get_column(i, j - 1)[k] + m_v_old.get_column(i + 1, j - 1)[k]);
          const double s_plus = X(i + 1, j + 1, k) + X(i, j + 1, k),
                       s_minus = X(i + 1, j - 1, k) + X(i, j - 1, k);
          U += mpdata::cross_term(dt, u, v_bar, s_plus, s_minus, s_plus + s_minus, dy, eps);

          if (dz_c > 0.0) {
            const double w_bar = 0.25 * (w_old(i, j, k) + w_old(i, j, k - 1) +
                                         w_old(i + 1, j, k) + w_old(i + 1, j, k - 1));
            const double c_plus = C(i, j, kp) + C(i + 1, j, kp),
                         c_minus = C(i, j, km) + C(i + 1, j, km);
            U += mpdata::cross_term(dt, u, w_bar, c_plus, c_minus, c_plus + c_minus, dz_c, eps);
          }
        }
        uf[k] = U;
      }

      // northern face
      {
        const double v = m_v_old.get_column(i, j)[k];
        double V = 0.0;
        if (v != 0.0) {
          V = mpdata::corrective_velocity(dt, v, X(i, j, k), X(i, j + 1, k), dy, eps);

          const double u_bar = 0.25 * (m_u_old.get_column(i, j + 1)[k] + m_u_old.get_column(i, j)[k] +
                                       m_u_old.get_column(i - 1, j)[k] + m_u_old.get_column(i - 1, j + 1)[k]);
          const double s_plus = X(i + 1, j + 1, k) + X(i + 1, j, k),
                       s_minus = X(i - 1, j + 1, k) + X(i - 1, j, k);
          V += mpdata::cross_term(dt, v, u_bar, s_plus, s_minus, s_plus + s_minus, dx, eps);

          if (dz_c > 0.0) {
            const double w_bar = 0.25 * (w_old(i, j, k) + w_old(i, j, k - 1) +
                                         w_old(i, j + 1, k) + w_old(i, j + 1, k - 1));
            const double c_plus = C(i, j, kp) + C(i, j + 1, kp),
                         c_minus = C(i, j, km) + C(i, j + 1, km);
            V += mpdata::cross_term(dt, v, w_bar, c_plus, c_minus, c_plus + c_minus, dz_c, eps);
          }
        }
        vf[k] = V;
      }

      // top face (between volumes k and k + 1)
      {
        double W = 0.0;
        if (k + 1 < Mz) {
          const double w = m_w_old.get_column(i, j)[k];
          if (w != 0.0) {
            W = mpdata::corrective_velocity(dt, w, C(i, j, k), C(i, j, k + 1), m_dz_face[k], eps);

            const double u_bar = 0.25 * (m_u_old.get_column(i, j)[k] + m_u_old.get_column(i - 1, j)[k] +
                                         m_u_old.get_column(i, j)[k + 1] + m_u_old.get_column(i - 1, j)[k + 1]);
            const double cx_plus = C(i + 1, j, k) + C(i + 1, j, k + 1),
                         cx_minus = C(i - 1, j, k) + C(i - 1, j, k + 1);
            W += mpdata::cross_term(dt, w, u_bar, cx_plus, cx_minus, cx_plus + cx_minus, dx, eps);

            const double v_bar = 0.25 * (m_v_old.get_column(i, j)[k] + m_v_old.get_column(i, j - 1)[k] +
                                         m_v_old.get_column(i, j)[k + 1] + m_v_old.get_column(i, j - 1)[k + 1]);
            const double cy_plus = C(i, j + 1, k) + C(i, j + 1, k + 1),
                         cy_minus = C(i, j - 1, k) + C(i, j - 1, k + 1);
            W += mpdata::cross_term(dt, w, v_bar, cy_plus, cy_minus, cy_plus + cy_minus, dy, eps);
          }
        }
        wf[k] = W;
      }
    }
  }
}

/*!
 * Flux-corrected-transport limiter of [Smolarkiewicz1990] applied to the anti-diffusive
 * velocities in `m_*_face`, using the result of the previous pass (`m_x_previous`) and
 * the input of the whole update (`m_x_input`) to define local bounds.
 */
void MPDATA3::limit(double dt, const array::Scalar1 &ice_thickness) {
  auto grid = m_x.grid();
  const int Mz = (int)grid->Mz();
  const double dx = grid->dx(), dy = grid->dy();

  // "up" and "down" limiting coefficients (beta) for every cell; ghosted so that the
  // coefficients of neighbors are available when scaling the face velocities
  array::Array3D beta_up(grid, "beta_up", array::WITH_GHOSTS, grid->z(), 1),
    beta_down(grid, "beta_down", array::WITH_GHOSTS, grid->z(), 1);

  {
    array::AccessScope list{ &ice_thickness, &m_x_previous, &m_x_input,
                             &m_u_face, &m_v_face, &m_w_face, &beta_up, &beta_down };

    auto X = [&](int i, int j, int k) -> double {
      return m_x_previous.get_column(i, j)[k];
    };
    auto X0 = [&](int i, int j, int k) -> double {
      return m_x_input.get_column(i, j)[k];
    };
    // concentration of the previous pass, in units of mass per unit area per meter
    auto C = [&](int i, int j, int k) -> double {
      double t = dz(k, ice_thickness(i, j));
      return t > 0.0 ? X(i, j, k) / t : 0.0;
    };

    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      double *b_up = beta_up.get_column(i, j), *b_down = beta_down.get_column(i, j);

      const double *uf = m_u_face.get_column(i, j), *uf_w = m_u_face.get_column(i - 1, j),
                   *vf = m_v_face.get_column(i, j), *vf_s = m_v_face.get_column(i, j - 1),
                   *wf = m_w_face.get_column(i, j);

      for (int k = 0; k < Mz; ++k) {
        const int kp = std::min(k + 1, Mz - 1), km = std::max(k - 1, 0);

        const double x_c = X(i, j, k), x_e = X(i + 1, j, k), x_w = X(i - 1, j, k),
                     x_n = X(i, j + 1, k), x_s = X(i, j - 1, k), x_t = X(i, j, kp),
                     x_b = X(i, j, km);

        double x_max = std::max({ x_c, x_e, x_w, x_n, x_s, x_t, x_b,
                                  X0(i, j, k), X0(i + 1, j, k), X0(i - 1, j, k),
                                  X0(i, j + 1, k), X0(i, j - 1, k), X0(i, j, kp), X0(i, j, km) });
        double x_min = std::min({ x_c, x_e, x_w, x_n, x_s, x_t, x_b,
                                  X0(i, j, k), X0(i + 1, j, k), X0(i - 1, j, k),
                                  X0(i, j + 1, k), X0(i, j - 1, k), X0(i, j, kp), X0(i, j, km) });

        const double u_e = uf[k], u_w = uf_w[k], v_n = vf[k], v_s = vf_s[k],
                     w_t = k + 1 < Mz ? wf[k] : 0.0, w_b = k > 0 ? wf[k - 1] : 0.0;

        // vertical fluxes transport the concentration
        const double c_c = C(i, j, k), c_t = k + 1 < Mz ? C(i, j, k + 1) : 0.0,
                     c_b = k > 0 ? C(i, j, k - 1) : 0.0;

        const double flux_in = dt * ((pp(u_w) * x_w - np(u_e) * x_e) / dx +
                                     (pp(v_s) * x_s - np(v_n) * x_n) / dy +
                                     (pp(w_b) * c_b - np(w_t) * c_t));
        const double flux_out = dt * ((pp(u_e) * x_c - np(u_w) * x_c) / dx +
                                      (pp(v_n) * x_c - np(v_s) * x_c) / dy +
                                      (pp(w_t) * c_c - np(w_b) * c_c));

        b_up[k]   = (x_max - x_c) / (flux_in + eps);
        b_down[k] = (x_c - x_min) / (flux_out + eps);
      }
    }
  }

  beta_up.update_ghosts();
  beta_down.update_ghosts();

  // scale the face velocities
  {
    array::AccessScope list{ &m_u_face, &m_v_face, &m_w_face, &beta_up, &beta_down };

    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      double *uf = m_u_face.get_column(i, j), *vf = m_v_face.get_column(i, j),
             *wf = m_w_face.get_column(i, j);

      const double *bu_c = beta_up.get_column(i, j), *bd_c = beta_down.get_column(i, j),
                   *bu_e = beta_up.get_column(i + 1, j), *bd_e = beta_down.get_column(i + 1, j),
                   *bu_n = beta_up.get_column(i, j + 1), *bd_n = beta_down.get_column(i, j + 1);

      for (int k = 0; k < Mz; ++k) {
        // east: flow from c to e limited by "down" at c and "up" at e, and vice versa
        if (uf[k] > 0.0) {
          uf[k] *= std::min(1.0, std::min(bd_c[k], bu_e[k]));
        } else if (uf[k] < 0.0) {
          uf[k] *= std::min(1.0, std::min(bu_c[k], bd_e[k]));
        }

        if (vf[k] > 0.0) {
          vf[k] *= std::min(1.0, std::min(bd_c[k], bu_n[k]));
        } else if (vf[k] < 0.0) {
          vf[k] *= std::min(1.0, std::min(bu_c[k], bd_n[k]));
        }

        if (k + 1 < Mz) {
          if (wf[k] > 0.0) {
            wf[k] *= std::min(1.0, std::min(bd_c[k], bu_c[k + 1]));
          } else if (wf[k] < 0.0) {
            wf[k] *= std::min(1.0, std::min(bu_c[k], bd_c[k + 1]));
          }
        }
      }
    }
  }
}

/*!
 * One explicit upwind step with the current face velocities: `m_x_previous` -> `m_x`.
 *
 * Outgoing fluxes of a cell are scaled down if they would remove more than the cell
 * contains during `dt`, which can happen in thin truncated top volumes even when the CFL
 * condition holds for full volumes. This keeps the solution non-negative and conserves
 * mass (both sides of a face use the same scaled flux).
 */
void MPDATA3::step(double dt, const array::Scalar1 &ice_thickness) {
  auto grid = m_x.grid();
  const int Mz = (int)grid->Mz();
  const double dx = grid->dx(), dy = grid->dy();

  array::Array3D scale(grid, "outflow_scale", array::WITH_GHOSTS, grid->z(), 1);

  // flux across the eastern face of (i, j) at level k, mass per unit area times velocity
  auto flux_x = [&](int i, int j, int k) -> double {
    return upwind(m_x_previous.get_column(i, j)[k], m_x_previous.get_column(i + 1, j)[k],
                  m_u_face.get_column(i, j)[k]);
  };
  auto flux_y = [&](int i, int j, int k) -> double {
    return upwind(m_x_previous.get_column(i, j)[k], m_x_previous.get_column(i, j + 1)[k],
                  m_v_face.get_column(i, j)[k]);
  };
  // flux across the top face of volume k (concentration times velocity)
  auto flux_z = [&](int i, int j, int k) -> double {
    if (k + 1 >= Mz) {
      return 0.0;
    }
    const double H = ice_thickness(i, j);
    const double t_c = dz(k, H), t_t = dz(k + 1, H);
    const double c_c = t_c > 0.0 ? m_x_previous.get_column(i, j)[k] / t_c : 0.0,
                 c_t = t_t > 0.0 ? m_x_previous.get_column(i, j)[k + 1] / t_t : 0.0;
    return upwind(c_c, c_t, m_w_face.get_column(i, j)[k]);
  };

  {
    array::AccessScope list{ &ice_thickness, &m_x_previous, &m_u_face, &m_v_face, &m_w_face, &scale };

    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      double *s = scale.get_column(i, j);

      for (int k = 0; k < Mz; ++k) {
        const double F_e = flux_x(i, j, k), F_w = flux_x(i - 1, j, k),
                     F_n = flux_y(i, j, k), F_s = flux_y(i, j - 1, k),
                     F_t = flux_z(i, j, k), F_b = k > 0 ? flux_z(i, j, k - 1) : 0.0;

        const double outflow = dt * ((pp(F_e) - np(F_w)) / dx + (pp(F_n) - np(F_s)) / dy +
                                     (pp(F_t) - np(F_b)));
        const double x_c = m_x_previous.get_column(i, j)[k];

        s[k] = (outflow > x_c and outflow > 0.0) ? std::max(x_c, 0.0) / outflow : 1.0;
      }
    }
  }

  scale.update_ghosts();

  {
    array::AccessScope list{ &ice_thickness, &m_x_previous, &m_u_face, &m_v_face, &m_w_face,
                             &scale, &m_x };

    // scaled flux: the scale factor of the cell the flux leaves
    auto limited = [&](double F, double s_from, double s_to) -> double {
      return F * (F >= 0.0 ? s_from : s_to);
    };

    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      const double *s_c = scale.get_column(i, j), *s_e = scale.get_column(i + 1, j),
                   *s_w = scale.get_column(i - 1, j), *s_n = scale.get_column(i, j + 1),
                   *s_s = scale.get_column(i, j - 1);
      const double *x_old = m_x_previous.get_column(i, j);
      double *x_new = m_x.get_column(i, j);

      for (int k = 0; k < Mz; ++k) {
        const double F_e = limited(flux_x(i, j, k), s_c[k], s_e[k]),
                     F_w = limited(flux_x(i - 1, j, k), s_w[k], s_c[k]),
                     F_n = limited(flux_y(i, j, k), s_c[k], s_n[k]),
                     F_s = limited(flux_y(i, j - 1, k), s_s[k], s_c[k]),
                     F_t = k + 1 < Mz ? limited(flux_z(i, j, k), s_c[k], s_c[k + 1]) : 0.0,
                     F_b = k > 0 ? limited(flux_z(i, j, k - 1), s_c[k - 1], s_c[k]) : 0.0;

        x_new[k] = x_old[k] - dt * ((F_e - F_w) / dx + (F_n - F_s) / dy + (F_t - F_b));
      }
    }
  }
}

void MPDATA3::update(double dt, const array::Scalar1 &ice_thickness,
                     const array::CellType1 &cell_type, const array::Array3D &x,
                     const array::Array3D &u, const array::Array3D &v,
                     const array::Array3D &w) {

  // ghosted copy of the input (needed by the limiter)
  m_x_input.copy_from(x);

  for (int n = 0; n < m_N; ++n) {
    if (n == 0) {
      m_x_previous.copy_from(x);
      compute_interface_velocity(ice_thickness, cell_type, u, v, w);
    } else {
      m_x_previous.copy_from(m_x);
      m_u_old.copy_from(m_u_face);
      m_v_old.copy_from(m_v_face);
      m_w_old.copy_from(m_w_face);
      compute_corrective_velocity(dt, ice_thickness);
      if (m_nonoscillatory) {
        // the limiter uses the (unlimited) anti-diffusive velocities in m_*_face and
        // scales them in place; it needs their ghosts
        m_u_face.update_ghosts();
        m_v_face.update_ghosts();
        m_w_face.update_ghosts();
        limit(dt, ice_thickness);
      }
    }

    m_u_face.update_ghosts();
    m_v_face.update_ghosts();
    m_w_face.update_ghosts();

    step(dt, ice_thickness);
  }
}

std::shared_ptr<TransportScheme3D> TransportScheme3D::create(std::shared_ptr<const Grid> grid,
                                                             const std::string &kind, int N,
                                                             bool nonoscillatory) {
  if (kind == "upwind") {
    return std::make_shared<MPDATA3>(grid, 1, false);
  }

  if (kind == "mpdata") {
    return std::make_shared<MPDATA3>(grid, N, nonoscillatory);
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "unknown 3D transport scheme '%s' (allowed: upwind, mpdata)",
                                kind.c_str());
}

} // end of namespace pism
