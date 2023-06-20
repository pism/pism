/* Copyright (C) 2023 PISM Authors
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
#include <cassert>
#include <gsl/gsl_interp.h>
#include <memory>

#include "pism/age/Isochrones.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Context.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"

/*!
 *
 * TO DO: merge consecutive layers if the sum of their maximum thicknesses is below a
 * threshold (possibly vertical grid resolution).
 *
 */

namespace pism {

Isochrones::Isochrones(std::shared_ptr<const Grid> grid)
  : Component(grid) {

  const auto& time = grid->ctx()->time();

  auto requested_times = m_config->get_string("isochrones.deposition_times");
  m_deposition_times = time->parse_times(requested_times);
  int n_max = m_config->get_number("isochrones.max_n_layers");

  if (m_deposition_times.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot process isochrones.deposition_times = '%s'",
                                  requested_times.c_str());
  }

  if (m_deposition_times.size() > (size_t)n_max) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "deposition times '%s' corresponds to %d isochronal layers, which exceeds the allowed maximum (%d)",
                                  requested_times.c_str(),
                                  (int)m_deposition_times.size(),
                                  n_max);
  }

  m_depths = std::make_shared<array::Array3D>(grid, "isochronal_layer_depths",
                                              array::WITHOUT_GHOSTS, m_deposition_times);

  m_depths->metadata().long_name("thicknesses of isochronal layers").units("m");
  auto &z = m_depths->metadata(0).z();
  z.clear_all_strings();
  z.set_name("deposition_time");
  z["units"] = time->units_string();
  z["calendar"] = time->calendar();
  z["long name"] = "minimum deposition time for an isochronal layer";
  z["axis"] = "T";

  m_tmp = std::make_shared<array::Array3D>(grid, "temporary storage",
                                           array::WITH_GHOSTS, m_deposition_times);
}

void Isochrones::init(const Geometry &geometry) {

  m_depths->set(0.0);

  // use ice thickness to create one layer
  //
  // FIXME: add several layers of equal thickness when bootstrapping (instead of just one)
  {
    array::AccessScope scope{&geometry.ice_thickness, m_depths.get()};

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      m_depths->get_column(i, j)[0] = geometry.ice_thickness(i, j);
    }

    m_top_layer = 0;
  }

  auto t = m_grid->ctx()->time()->current();

  if (t < m_deposition_times[0]) {
    m_time_index = -1;
  } else {
    int N = static_cast<int>(m_deposition_times.size());

    // Find the index k such that m_deposition_times[k] <= T < m_deposition_times[k + 1]
    //
    // Note: `m_time_index` below will be strictly less than `N - 1`.
    m_time_index = gsl_interp_bsearch(m_deposition_times.data(), t, 0, N - 1);
  }

  // read from an input file when restarting
  //
  // FIXME: implement restarting
}

void Isochrones::update(double t, double dt,
                        const array::Array3D &u,
                        const array::Array3D &v,
                        const array::Scalar &ice_thickness,
                        const array::Scalar &climatic_mass_balance,
                        const array::Scalar &basal_melt_rate) {

  auto ice_density = m_config->get_number("constants.ice.density");

  {
    array::AccessScope scope{&climatic_mass_balance, &basal_melt_rate, m_depths.get()};

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double *d = m_depths->get_column(i, j);

      // apply the surface mass balance
      {
        double dH = dt * climatic_mass_balance(i, j) / ice_density;

        // apply thickness change to a layer, starting from the top-most
        for (int k = m_top_layer; k >= 0; --k) {
          if (d[k] + dH >= 0.0) {
            // thickness change is non-negative or does not remove the whole layer: apply to
            // the current layer and stop
            d[k] += dH;
            break;
          }

          dH += d[k];
          d[k] = 0.0;
        }
      }
      // apply the basal melt rate
      {
        double dH = -dt * basal_melt_rate(i, j) / ice_density;

        // apply thickness change to a layer, starting from the bottom
        for (int k = 0; k <= m_top_layer; ++k) {
          if (d[k] + dH >= 0.0) {
            // thickness change is non-negative or does not remove the whole layer: apply to
            // the current layer and stop
            d[k] += dH;
            break;
          }

          dH += d[k];
          d[k] = 0.0;
        }
      }
    }
  }

  // note: this updates ghosts of m_tmp
  m_tmp->copy_from(*m_depths);

  array::AccessScope scope{&u, &v, m_depths.get(), m_tmp.get(), &ice_thickness};

  double
    dx = m_grid->dx(),
    dy = m_grid->dy(),
    H_min = m_config->get_number("geometry.ice_free_thickness_standard");

  // flux estimated using first-order upwinding
  auto Q = [](double U, double f_n, double f_p) {
    return U * (U >= 0 ? f_n : f_p);
  };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      *d_c   = m_tmp->get_column(i, j),
      *d_n = m_tmp->get_column(i, j + 1),
      *d_e = m_tmp->get_column(i + 1, j),
      *d_s = m_tmp->get_column(i, j - 1),
      *d_w = m_tmp->get_column(i - 1, j);

    double *d = m_depths->get_column(i, j);

    pism::stencils::Star<double> z = 0.0;
    double d_total = 0.0;
    for (int k = 0; k <= m_top_layer; ++k) {

      // Evaluate velocities in the *middle* (vertically) of the current layer. I am
      // guessing that in many applications near the base of the ice layers get thin, so
      // *sampling* is okay because a layer thickness is likely to be smaller than the
      // vertical grid resolution used by the stress balance model. In the upper half of
      // ice thickness, on the other hand, we have less variation of ice velocity in the
      // vertical (du/dz is smaller), so we don't lose too much accuracy by linearizing
      // u(z). For a linear f(x) we have
      //
      // (f(a) + f(b))/2 = f((a + b) / 2) = f(a + (b - a) / 2),
      //
      // which allows us to estimate the "average horizontal velocity in a layer" using
      // *one* interpolation.
      //
      // Note, however, that the modeled depth of an isochrone is affected by depths of
      // all the isochrones below it since the elevation of an isochrone above the bed is
      // the sum of depths of all the layers below it.
      //
      // This implies that we should have at least a few layers *below* an isochrone we're
      // interested in.

      double
        U   = u.interpolate(i, j, z.c + 0.5 * d_c[k]),
        U_e = u.interpolate(i + 1, j, z.e + 0.5 * d_e[k]),
        U_w = u.interpolate(i - 1, j, z.w + 0.5 * d_w[k]);

      double
        V   = v.interpolate(i, j, z.c + 0.5 * d_c[k]),
        V_n = v.interpolate(i, j + 1, z.n + 0.5 * d_n[k]),
        V_s = v.interpolate(i, j - 1, z.s + 0.5 * d_s[k]);

      double
        Q_n = Q(0.5 * (V + V_n), d_c[k], d_n[k]),
        Q_e = Q(0.5 * (U + U_e), d_c[k], d_e[k]),
        Q_s = Q(0.5 * (V + V_s), d_s[k], d_c[k]),
        Q_w = Q(0.5 * (U + U_w), d_w[k], d_c[k]);

      d[k] = d_c[k] - dt * ((Q_e - Q_w) / dx + (Q_n - Q_s) / dy);

      assert(d[k] >= 0.0);

      // ensure non-negativity (should not be necessary, but still)
      d[k] = std::max(d[k], 0.0);

      d_total += d[k];

      z.c += d_c[k];
      z.n += d_n[k];
      z.e += d_e[k];
      z.s += d_s[k];
      z.w += d_w[k];
    }

    assert(ice_thickness(i, j) < H_min or d_total > 0.0);

    // re-scale so that the sum of layer thicknesses is equal to the ice_thickness
    if (d_total > 0.0) {
      double S = ice_thickness(i, j) / d_total;
      for (int k = 0; k <= m_top_layer; ++k) {
        d[k] *= S;
      }
    }
  }

  // add one more layer if we reached the next deposition time
  {
    double T = t + dt;
    int N = m_deposition_times.size();

    // Find the index k such that m_deposition_times[k] <= T < m_deposition_times[k + 1]
    //
    // Note: `k` below will be strictly less than `N - 1`.
    //
    // FIXME: consider using a gsl_interp_accel to speed this up
    size_t k = gsl_interp_bsearch(m_deposition_times.data(), T, 0, N - 1);

    if (k > (size_t)m_time_index and m_top_layer < N - 1) {

      m_top_layer += 1;
      m_time_index = k;

      if (m_top_layer == N - 1) {
        m_log->message(2, "reached isochrones.max_n_layers\n");
      }
    }
  }
}

MaxTimestep Isochrones::max_timestep_impl(double t) const {
  double t0 = m_deposition_times[0];
  if (t < t0) {
    return {t0 - t, "isochrones"};
  }

  if (t >= m_deposition_times.back()) {
    return {"isochrones"};
  }

  auto N = m_deposition_times.size();

  // Find the index k such that m_deposition_times[k] <= T < m_deposition_times[k + 1]
  //
  // Note: `k` below will be strictly less than `N - 1`.
  //
  // FIXME: use a gsl_interp_accel to speed this up
  size_t k = gsl_interp_bsearch(m_deposition_times.data(), t, 0, N - 1);

  return {m_deposition_times[k + 1] - t, "isochrones"};
}

void Isochrones::define_model_state_impl(const File &output) const {
  m_depths->define(output, io::PISM_DOUBLE);
}

void Isochrones::write_model_state_impl(const File &output) const {
  m_depths->write(output);
}

const array::Array3D& Isochrones::layer_depths() const {
  return *m_depths;
}

namespace diagnostics {

/*! @brief Report isochrone depth */
class IsochroneDepths : public Diag<Isochrones>
{
public:
  IsochroneDepths(const Isochrones *m)
    : Diag<Isochrones>(m) {

    const auto& time = m_grid->ctx()->time();

    m_vars = {{m_sys, "isochrone_depth", model->layer_depths().levels()}};

    m_vars[0].long_name("isochrone depth").units("m");
    auto &z = m_vars[0].z();
    z.clear_all_strings();
    z.set_name("deposition_time");
    z["units"] = time->units_string();
    z["calendar"] = time->calendar();
    z["long name"] = "minimum deposition time for an isochronal layer";
    z["axis"] = "T";
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {

    const auto& layer_depths = model->layer_depths();

    auto result = layer_depths.duplicate();
    result->metadata(0) = m_vars[0];

    int N = result->levels().size();

    array::AccessScope scope{&layer_depths, result.get()};

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double *column = result->get_column(i, j);
      const double *d = layer_depths.get_column(i, j);

      double total_depth = 0.0;
      for (int k = N - 1; k >= 0; --k) {
        total_depth += d[k];
        column[k] = total_depth;
      }
    }

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList Isochrones::diagnostics_impl() const {
  return {{"isochrone_depth", Diagnostic::Ptr(new diagnostics::IsochroneDepths(this))}};
}

} // end of namespace pism
