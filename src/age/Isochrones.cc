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
#include <alloca.h>
#include <cassert>
#include <cstddef>
#include <gsl/gsl_interp.h>
#include <memory>
#include <vector>

#include "pism/age/Isochrones.hh"
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

namespace details {
static const char *layer_count_variable_name     = "isochronal_layer_count";
static const char *layer_thickness_variable_name = "isochronal_layer_thickness";
static const char *deposition_time_variable_name = "deposition_time";
static const char *isochrone_depth_variable_name = "isochrone_depth";

static const char *times_parameter = "isochrones.deposition_times";
static const char *N_max_parameter = "isochrones.max_n_layers";

static const char *N_boot_parameter = "isochrones.bootstrapping.n_layers";
}

//
/*!
 * Allocates storage and initializes
 *
 * - requested deposition times
 *
 * but does not initialize
 *
 * - layer thicknesses
 * - deposition times for existing layers
 * - topmost layer index
 */
Isochrones::Isochrones(std::shared_ptr<const Grid> grid)
  : Component(grid) {

  using namespace details;

  auto N_max = (int)m_config->get_number(N_max_parameter);

  if (N_max < 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s have to be non-negative (got %d)",
                                  N_max_parameter, N_max);
  }

  const auto &time = grid->ctx()->time();

  auto requested_times = m_config->get_string(times_parameter);
  m_deposition_times   = time->parse_times(requested_times);

  auto N_deposition_times = m_deposition_times.size();
  if (N_deposition_times == 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot process isochrones.deposition_times = '%s'",
                                  requested_times.c_str());
  }

  if ((int)N_deposition_times > N_max) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION,
        "the number of times %d in '%s' exceeds the amount of storage allocated (%s = %d)",
        (int)N_deposition_times, times_parameter, N_max_parameter, (int)N_max);
  }

  // Note: array::Array delays allocation until the last moment, so we can cheaply
  // re-allocate storage if the number of "levels" used here turns out to be
  // inappropriate.
  allocate(m_deposition_times);
}


void Isochrones::allocate(const std::vector<double> &levels) {
  using namespace details;

  const auto &time = m_grid->ctx()->time();

  m_layer_thickness = std::make_shared<array::Array3D>(m_grid, layer_thickness_variable_name,
                                                       array::WITHOUT_GHOSTS, levels);

  m_layer_thickness->metadata().long_name("thicknesses of isochronal layers").units("m");

  auto z_description =
      pism::printf("times for isochrones in '%s'; earliest deposition times for layers in '%s'",
                   isochrone_depth_variable_name, layer_thickness_variable_name);
  auto &z = m_layer_thickness->metadata(0).z();
  z.clear()
      .set_name(deposition_time_variable_name)
      .long_name(z_description)
      .units(time->units_string());
  z["calendar"] = time->calendar();

  m_tmp = std::make_shared<array::Array3D>(m_grid, "temporary storage", array::WITH_GHOSTS, levels);
}

/*!
 * When bootstrapping, we can put all the existing ice thickness into the bottom layer and
 * then keep adding to it until we reach the next deposition time, *or*, we can distribute
 * the ice thickness among N "bootstrapping" layers and then apply SMB to the layer `N+1`.
 * The second option allows us to increase accuracy: the quality of the horizontal
 * velocity approximation used to transport mass within layers is higher if the layers are
 * thin.
 */
void Isochrones::bootstrap(const array::Scalar &ice_thickness) {
  using namespace details;

  m_layer_thickness->set(0.0);

  auto N_bootstrap        = static_cast<int>(m_config->get_number(N_boot_parameter));
  auto N_max              = static_cast<int>(m_config->get_number(N_max_parameter));
  auto N_deposition_times = m_deposition_times.size();

  if (N_bootstrap < 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s have to be non-negative (got %d)",
                                  N_boot_parameter, N_bootstrap);
  }

  m_log->message(2, "* Bootstrapping the isochrone tracking model, adding %d isochronal layers...\n",
                 N_bootstrap);

  if (N_bootstrap + (int)N_deposition_times > N_max) {
    auto deposition_times = m_config->get_string("isochrones.deposition_times");
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION, "%s (%d) + %s (%d) exceeds the amount of storage allocated (%s = %d)",
        N_boot_parameter, (int)N_bootstrap, times_parameter,
        (int)N_deposition_times, N_max_parameter, (int)N_max);
  }

  m_top_layer_index = 0;
  if (N_bootstrap > 0) {
    // prepend "bootstrapping" layers to m_deposition_times
    std::vector<double> requested_deposition_times = m_deposition_times;
    double T_0                                     = m_grid->ctx()->time()->current();

    m_deposition_times.clear();
    for (int k = 0; k < N_bootstrap; ++k) {
      m_deposition_times.push_back(T_0);
      m_top_layer_index += 1;
    }
    for (const auto &t : requested_deposition_times) {
      m_deposition_times.push_back(t);
    }

    // re-allocate storage
    allocate(m_deposition_times);

    array::AccessScope scope{ &ice_thickness, m_layer_thickness.get() };

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double H = ice_thickness(i, j);

      double *column = m_layer_thickness->get_column(i, j);
      for (int k = 0; k < N_bootstrap; ++k) {
        column[k] = H / static_cast<double>(N_bootstrap);
      }
    }
  } else {
    array::AccessScope scope{ &ice_thickness, m_layer_thickness.get() };

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();
      m_layer_thickness->get_column(i, j)[0] = ice_thickness(i, j);
    }
  }
}

void Isochrones::restart(const File &input_file, int record) {
  using namespace details;

  m_log->message(2, "* Initializing the isochrone tracking model from '%s'...\n",
                 input_file.filename().c_str());


  // get deposition times from the input file
  std::vector<double> old_deposition_times;
  {
    auto n_deposition_times = input_file.dimension_length(deposition_time_variable_name);

    // read
    old_deposition_times.resize(n_deposition_times);
    input_file.read_variable(deposition_time_variable_name, { 0 }, { n_deposition_times },
                             old_deposition_times.data());
  }

  // Add requested deposition times
  //
  // This trickery is needed because "-isochrones.deposition_times 1000" will generate
  // deposition times every 1000 years for the duration of the current run... and when we
  // are re-starting we need to include times from both the current and the /previous/ run.
  double last_time                         = old_deposition_times.back();
  std::vector<double> new_deposition_times = old_deposition_times;
  for (auto t : m_deposition_times) {
    if (t > last_time) {
      new_deposition_times.push_back(t);
    }
  }

  // check if we are allowed to allocate storage for this many layers
  auto N_max = (size_t)m_config->get_number(N_max_parameter);

  if (new_deposition_times.size() > N_max) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION,
        "the total number of isochronal layers (from the input file '%s' plus requested) exceeds '%s' = %d",
        input_file.filename().c_str(), N_max_parameter, (int)N_max);
  }

  m_deposition_times = new_deposition_times;

  // re-allocate storage
  allocate(m_deposition_times);

  // allocate temporary storage, read in layer thicknesses, move layer thicknesses from
  // temporary storage into m_layer_thickness:
  {
    auto tmp = std::make_shared<array::Array3D>(m_grid, layer_thickness_variable_name,
                                                array::WITHOUT_GHOSTS, old_deposition_times);
    tmp->metadata().long_name("thicknesses of isochronal layers").units("m");
    auto &z = tmp->metadata().z();
    z.clear()
        .set_name(deposition_time_variable_name)
        .units(m_layer_thickness->metadata().z()["units"]);

    tmp->read(input_file, record);

    array::AccessScope scope{ tmp.get(), m_layer_thickness.get() };

    size_t N = tmp->get_levels().size();
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto *input  = tmp->get_column(i, j);
      auto *output = m_layer_thickness->get_column(i, j);

      for (size_t k = 0; k < N; ++k) {
        output[k] = input[k];
      }
    }
  }

  // set m_top_layer_index
  {
    double n_active_layers = 0;
    input_file.read_variable(layer_count_variable_name, { (unsigned int)record }, { 1 },
                             &n_active_layers);

    m_top_layer_index = static_cast<size_t>(n_active_layers) - 1;
  }
}

void Isochrones::update(double t, double dt, const array::Array3D &u, const array::Array3D &v,
                        const array::Scalar &ice_thickness,
                        const array::Scalar &climatic_mass_balance,
                        const array::Scalar &basal_melt_rate) {

  auto ice_density = m_config->get_number("constants.ice.density");

  // apply top surface and basal mass balance terms:
  {
    array::AccessScope scope{ &climatic_mass_balance, &basal_melt_rate, m_layer_thickness.get() };

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double *d = m_layer_thickness->get_column(i, j);

      // apply the surface mass balance
      {
        double dH = dt * climatic_mass_balance(i, j) / ice_density;

        // apply thickness change to a layer, starting from the top-most
        for (int k = (int)m_top_layer_index; k >= 0; --k) {
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
        for (size_t k = 0; k <= m_top_layer_index; ++k) {
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

  // transport mass within layers:
  {
    // note: this updates ghosts of m_tmp
    m_tmp->copy_from(*m_layer_thickness);

    array::AccessScope scope{ &u, &v, m_layer_thickness.get(), m_tmp.get(), &ice_thickness };

    double dx = m_grid->dx(), dy = m_grid->dy(),
           H_min = m_config->get_number("geometry.ice_free_thickness_standard");

    // flux estimated using first-order upwinding
    auto Q = [](double U, double f_n, double f_p) { return U * (U >= 0 ? f_n : f_p); };

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *d_c = m_tmp->get_column(i, j), *d_n = m_tmp->get_column(i, j + 1),
                   *d_e = m_tmp->get_column(i + 1, j), *d_s = m_tmp->get_column(i, j - 1),
                   *d_w = m_tmp->get_column(i - 1, j);

      double *d = m_layer_thickness->get_column(i, j);

      pism::stencils::Star<double> z = 0.0;
      double d_total                 = 0.0;
      for (size_t k = 0; k <= m_top_layer_index; ++k) {

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

        double U   = u.interpolate(i, j, z.c + 0.5 * d_c[k]),
               U_e = u.interpolate(i + 1, j, z.e + 0.5 * d_e[k]),
               U_w = u.interpolate(i - 1, j, z.w + 0.5 * d_w[k]);

        double V   = v.interpolate(i, j, z.c + 0.5 * d_c[k]),
               V_n = v.interpolate(i, j + 1, z.n + 0.5 * d_n[k]),
               V_s = v.interpolate(i, j - 1, z.s + 0.5 * d_s[k]);

        double Q_n = Q(0.5 * (V + V_n), d_c[k], d_n[k]), Q_e = Q(0.5 * (U + U_e), d_c[k], d_e[k]),
               Q_s = Q(0.5 * (V + V_s), d_s[k], d_c[k]), Q_w = Q(0.5 * (U + U_w), d_w[k], d_c[k]);

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
        for (size_t k = 0; k <= m_top_layer_index; ++k) {
          d[k] *= S;
        }
      }
    }
  }
  // add one more layer if we reached the next deposition time
  {
    double T                  = t + dt;
    size_t N_deposition_times = m_deposition_times.size();
    size_t max_n_levels       = m_layer_thickness->get_levels().size();

    // Find the index k such that m_deposition_times[k] <= T < m_deposition_times[k + 1]
    //
    // Note: `k` below will be strictly less than `N - 1`, ensuring that the index "k + 1"
    // is valid.
    //
    // FIXME: consider using a gsl_interp_accel to speed this up
    size_t k = gsl_interp_bsearch(m_deposition_times.data(), T, 0, N_deposition_times - 1);

    double T_k = m_deposition_times[k];

    double top_layer_deposition_time = m_layer_thickness->get_levels().at(m_top_layer_index);
    if (T_k > top_layer_deposition_time) {
      // we reached the next requested deposition time

      if (m_top_layer_index < max_n_levels - 1) {
        // not too many layers yet: add one more
        m_top_layer_index += 1;
      } else {
        // we have as many layers as we can handle: keep adding to the top layer
        m_log->message(2,
                       "Isochrone tracking: reached isochrones.max_n_layers and can't add more.\n"
                       "  SMB will contribute to the current top layer.");
      }
    }
  }
}


MaxTimestep Isochrones::max_timestep_impl(double t) const {
  double t0 = m_deposition_times[0];
  if (t < t0) {
    return { t0 - t, "isochrones" };
  }

  if (t >= m_deposition_times.back()) {
    return { "isochrones" };
  }

  auto N = m_deposition_times.size();

  // Find the index k such that m_deposition_times[k] <= T < m_deposition_times[k + 1]
  //
  // Note: `k` below will be strictly less than `N - 1`, ensuring that the index "k + 1"
  // is valid.
  //
  // FIXME: consider using gsl_interp_accel to speed this up
  size_t k = gsl_interp_bsearch(m_deposition_times.data(), t, 0, N - 1);

  return { m_deposition_times[k + 1] - t, "isochrones" };
}

void Isochrones::define_model_state_impl(const File &output) const {
  using namespace details;

  m_layer_thickness->define(output, io::PISM_DOUBLE);

  if (not output.find_variable(layer_count_variable_name)) {
    auto time_name = m_config->get_string("time.dimension_name");

    output.define_variable(layer_count_variable_name, io::PISM_INT, { time_name });
    output.write_attribute(layer_count_variable_name, "long_name",
                           pism::printf("number of 'active' isochronal layers in '%s'",
                                        m_layer_thickness->metadata().get_name().c_str()));
  }
}

void Isochrones::write_model_state_impl(const File &output) const {
  m_layer_thickness->write(output);

  double N = static_cast<double>(m_top_layer_index) + 1;

  unsigned int t_last = output.nrecords() - 1;

  output.write_variable(details::layer_count_variable_name, { t_last }, { 1 }, &N);
}

const array::Array3D &Isochrones::layer_thicknesses() const {
  return *m_layer_thickness;
}

namespace diagnostics {

/*! @brief Report isochrone depth */
class IsochroneDepths : public Diag<Isochrones> {
public:
  IsochroneDepths(const Isochrones *m) : Diag<Isochrones>(m) {
    using namespace details;

    const auto &time = m_grid->ctx()->time();

    m_vars = { { m_sys, isochrone_depth_variable_name, model->layer_thicknesses().get_levels() } };

    auto description = pism::printf("depth below surface of isochrones for times in '%s'",
                                    deposition_time_variable_name);

    m_vars[0].long_name(description).units("m");
    auto &z = m_vars[0].z();
    z.clear()
        .set_name(deposition_time_variable_name)
        .long_name(
            pism::printf("deposition times for isochrones in '%s'", isochrone_depth_variable_name))
        .units(time->units_string());
    z["calendar"] = time->calendar();
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {

    const auto &layer_thicknesses = model->layer_thicknesses();

    auto result         = layer_thicknesses.duplicate();
    result->metadata(0) = m_vars[0];

    size_t N = result->get_levels().size();

    array::AccessScope scope{ &layer_thicknesses, result.get() };

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double *column  = result->get_column(i, j);
      const double *d = layer_thicknesses.get_column(i, j);

      double total_depth = 0.0;
      for (int k = (int)N - 1; k >= 0; --k) {
        total_depth += d[k];
        column[k] = total_depth;
      }
    }

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList Isochrones::diagnostics_impl() const {
  return { { details::isochrone_depth_variable_name,
             Diagnostic::Ptr(new diagnostics::IsochroneDepths(this)) } };
}

} // end of namespace pism

/*
 * LocalWords: LocalWords deposition
*/
