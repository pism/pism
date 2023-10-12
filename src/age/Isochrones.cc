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
#include <cassert>              // assert
#include <cstddef>              // size_t
#include <gsl/gsl_interp.h>
#include <memory>               // std::shared_ptr
#include <vector>

#include "pism/age/Isochrones.hh"
#include "pism/util/Context.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/interpolation.hh"
#include "pism/util/petscwrappers/Vec.hh"

/*!
 * Initialization cases:
 *
 * Parse requested deposition times and discard times outside of the modeled time span.
 * Store these in `times`.
 *
 * 1. Restarting.
 *
 * - Read deposition times `input_times` from the input file. Find all the `input_times`
 *   that are *before the beginning of the modeled time span* and record their indexes in
 *   `input_times`; save these to `old_times`.
 *
 * - Append `times` to `old_times`. It's okay if `times` is empty.
 *
 * - Read layer thicknesses corresponding to all `old_times` from the input file using
 *   recorded indexes.
 *
 * - Allocate storage for layer thicknesses for combined deposition times.
 *
 * - Renormalize layer thicknesses to ensure that their sum adds up to the ice thickness.
 *
 * 2. Restarting with regridding.
 *
 * - Same as restarting, but using the regridding file instead of the input file and
 *   linear interpolation to read in layer thicknesses.
 *
 * 3. Bootstrapping.
 *
 * - If `isochrones.bootstrapping.n_layers` is zero, create one layer with the thickness
 *   equal to the ice thickness. The SMB will be applied to this layer until the model
 *   reaches a deposition time and a new layer is created.
 *
 * - If `N = isochrones.bootstrapping.n_layers` is positive, create `N` layers with
 *   thicknesses equal to the ice thickness divided by `N`. Then add one more layer of
 *   zero thickness. The SMB will be applied to this layer until the model reaches a
 *   deposition time and a new layer is created.
 *
 * - Set deposition times to `times`. Stop if `times` is empty.
 *
 * 4. Bootstrapping with regridding.
 *
 * - Same as restarting with regridding.
 *
 */

namespace pism {

namespace details {
static const char *layer_thickness_variable_name = "isochronal_layer_thickness";
static const char *deposition_time_variable_name = "deposition_time";
static const char *isochrone_depth_variable_name = "isochrone_depth";

static const char *times_parameter = "isochrones.deposition_times";
static const char *N_max_parameter = "isochrones.max_n_layers";
static const char *N_boot_parameter = "isochrones.bootstrapping.n_layers";


//! Checks if a vector of doubles is not decreasing.
static bool non_decreasing(const std::vector<double> &a) {
  size_t N = a.size();

  if (N < 2) {
    return true;
  }

  for (size_t k = 0; k < N - 1; k++) {
    if (a[k] > a[k + 1]) {
      return false;
    }
  }

  return true;
}

/*!
 * Allocate storage for layer thicknesses given a list of deposition times `times`.
 *
 * Also sets metadata.
 */
static std::shared_ptr<array::Array3D> allocate_layer_thickness(std::shared_ptr<const Grid> grid,
                                                                const std::vector<double> &times) {
  using namespace details;

  auto N_max = (int)grid->ctx()->config()->get_number(N_max_parameter);
  if ((int)times.size() > N_max) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "the total number of isochronal layers (%d) exceeds '%s' = %d",
                                  (int)times.size(), N_max_parameter, (int)N_max);
  }

  const auto &time = grid->ctx()->time();

  auto result = std::make_shared<array::Array3D>(grid, layer_thickness_variable_name,
                                                 array::WITHOUT_GHOSTS, times);

  result->metadata().long_name("thicknesses of isochronal layers").units("m");

  auto z_description =
      pism::printf("times for isochrones in '%s'; earliest deposition times for layers in '%s'",
                   isochrone_depth_variable_name, layer_thickness_variable_name);
  auto &z = result->metadata(0).z();
  z.clear()
      .set_name(deposition_time_variable_name)
      .long_name(z_description)
      .units(time->units_string());
  z["calendar"] = time->calendar();

  return result;
}

/*!
 * Allocate layer thicknesses and copy relevant info from `input`.
 *
 * @param[in] input input layer thicknesses and deposition times, read from an input file
 * @param[in] T_start start time of the current run
 * @param[in] requested_times requested deposition times
 */
static std::shared_ptr<array::Array3D>
allocate_layer_thickness(const array::Array3D &input, double T_start,
                         const std::vector<double> &requested_times) {

  auto grid = input.grid();

  const auto &input_times = input.get_levels();

  // the sequence of deposition times has to be monotonically non-decreasing
  // (bootstrapping may create several layers with identical deposition times)
  if (!details::non_decreasing(input_times)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "deposition times in '%s' have to be non-decreasing",
                                  input.get_name().c_str());
  }

  std::vector<double> deposition_times;
  for (auto t : input_times) {
    if (t <= T_start) {
      deposition_times.push_back(t);
    }
  }
  auto N_layers_to_copy = deposition_times.size();

  double last_kept_time = deposition_times.back();
  for (auto t : requested_times) {
    if (t > last_kept_time) {
      deposition_times.push_back(t);
    }
  }

  auto result = allocate_layer_thickness(grid, deposition_times);

  array::AccessScope scope{ &input, result.get() };

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const auto *in = input.get_column(i, j);
    auto *out      = result->get_column(i, j);

    for (size_t k = 0; k < N_layers_to_copy; ++k) {
      out[k] = in[k];
    }
  }

  return result;
}

/*!
 * Discard requested deposition times before the beginning and after the end of the run.
 */
static std::vector<double> prune_deposition_times(const Time &time,
                                                  const std::vector<double> &times) {
  double T_start = time.start(), T_end = time.end();

  std::vector<double> result;
  for (auto t : times) {
    if (t >= T_start and t <= T_end) {
      result.push_back(t);
    }
  }


  return result;
}

/*!
 * Read deposition times from an `input_file`.
 */
static std::vector<double> deposition_times(const File &input_file) {
  using namespace details;

  auto n_deposition_times = input_file.dimension_length(deposition_time_variable_name);

  std::vector<double> result(n_deposition_times);
  input_file.read_variable(deposition_time_variable_name, { 0 }, { n_deposition_times },
                           result.data());

  return result;
}

/*!
 * Process configuration parameters and return requested deposition times, excluding times
 * outside the current modeled time interval.
 */
std::vector<double> deposition_times(const Config &config, const Time &time) {
  using namespace details;
  try {
    auto N_max = (int)config.get_number(N_max_parameter);

    if (N_max < 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s have to be non-negative (got %d)",
                                    N_max_parameter, N_max);
    }

    auto requested_times = config.get_string(times_parameter);

    auto deposition_times = prune_deposition_times(time, time.parse_times(requested_times));

    auto N_deposition_times = deposition_times.size();
    if (N_deposition_times == 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "'%s' = '%s' has times within the modeled time span [%s, %s]",
                                    times_parameter, requested_times.c_str(),
                                    time.date(time.start()).c_str(), time.date(time.end()).c_str());
    }

    if ((int)N_deposition_times > N_max) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "the number of times (%d) in '%s' exceeds the amount of storage allocated ('%s' = %d)",
          (int)N_deposition_times, times_parameter, N_max_parameter, N_max);
    }

    return deposition_times;
  } catch (RuntimeError &e) {
    e.add_context("processing parameter '%s'", times_parameter);
    throw;
  }
}

/*!
 * Read layer deposition times and layer thicknesses from `input_file`.
 *
 * Uses bilinear interpolation in the X and Y directions.
 */
static std::shared_ptr<array::Array3D> regrid_layer_thickness(std::shared_ptr<const Grid> grid,
                                                              const File &input_file, int record) {

  auto result = details::allocate_layer_thickness(grid, deposition_times(input_file));

  auto N = (int)result->get_levels().size();

  auto metadata = result->metadata(0);

  auto variable_info = input_file.find_variable(metadata.get_name(), metadata["standard_name"]);

  grid::InputGridInfo input_grid(input_file, variable_info.name, metadata.unit_system(),
                                 grid->registration());

  // Set up 2D interpolation:
  LocalInterpCtx lic(input_grid, *grid, LINEAR);
  lic.start[T_AXIS] = record;
  lic.count[T_AXIS] = 1;
  lic.start[Z_AXIS] = 0;
  lic.count[Z_AXIS] = N;
  // Create an "identity" version of the interpolation in the "vertical" direction:
  {
    std::vector<double> Z(N);
    for (int k = 0; k < N; ++k) {
      Z[k] = k;
    }

    // note matching input and output grids below:
    lic.z = std::make_shared<Interpolation>(NEAREST, Z, Z);
  }

  petsc::VecArray tmp(result->vec());
  io::regrid_spatial_variable(metadata, *grid, lic, input_file, tmp.get());

  return result;
}

/*!
 * Read layer thickness from an `input_file` (without interpolation).xo
 */
static std::shared_ptr<array::Array3D> read_layer_thickness(std::shared_ptr<const Grid> grid,
                                                            const File &input_file, int record) {
  auto result = allocate_layer_thickness(grid, deposition_times(input_file));

  result->read(input_file, record);

  return result;
}

/*!
 * Compute the number of "active" layers in `deposition_times`.
 *
 * A layer is "active" if the model reached its deposition time.
 */
static size_t n_active_layers(std::vector<double> deposition_times, double current_time) {
  size_t result = 0;
  for (auto t : deposition_times) {
    if (t <= current_time) {
      result += 1;
    }
  }

  return result;
}

/*!
 * Returns `true` if the user asked to regrid layer thicknesses, otherwise `false`.
 */
static bool regridp(const Config &config) {

  auto regrid_file = config.get_string("input.regrid.file");

  if (regrid_file.empty()) {
    return false;
  }

  auto regrid_vars = set_split(config.get_string("input.regrid.vars"), ',');

  return member(details::layer_thickness_variable_name, regrid_vars);
}

/*!
 * Re-scale layer thicknesses so that the sum of layer thicknesses is equal to the
 * ice_thickness.
 *
 * @param[in] ice_thickness ice thickness, meters
 * @param[in,out] layer_thickness isochronal layer thickness
 */
static void renormalize(const array::Scalar &ice_thickness,
                        array::Array3D &layer_thickness) {
  auto grid = layer_thickness.grid();

  auto N = layer_thickness.get_levels().size();

  array::AccessScope scope{&ice_thickness, &layer_thickness};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *H = layer_thickness.get_column(i, j);

    double H_total = 0.0;
    for (int k = 0; k < (int)N; ++k) {
      H_total += H[k];
    }

    // re-scale so that the sum of layer thicknesses is equal to the ice_thickness
    if (H_total > 0.0) {
      double S = ice_thickness(i, j) / H_total;
      for (size_t k = 0; k < N; ++k) {
        H[k] *= S;
      }
    }
  }
}

} // namespace details


Isochrones::Isochrones(std::shared_ptr<const Grid> grid,
                       std::shared_ptr<const stressbalance::StressBalance> stress_balance)
    : Component(grid), m_stress_balance(stress_balance) {

  auto time = grid->ctx()->time();

  // Allocate storage for *one* isochronal layer just to ensure that this instance is in
  // a valid state when the constrictor returns.
  //
  // Note: array::Array delays allocation until the last moment, so we can cheaply
  // re-allocate storage if the number of "levels" used here turns out to be
  // inappropriate.
  m_layer_thickness = details::allocate_layer_thickness(m_grid, { time->current() });
  m_tmp             = m_layer_thickness->duplicate(array::WITH_GHOSTS);
  m_top_layer_index = 0;
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
  try {
    if (regridp(*m_config)) {
      File regrid_file(m_grid->com, m_config->get_string("input.regrid.file"), io::PISM_GUESS,
                       io::PISM_READONLY);

      auto last_record =
          regrid_file.nrecords(details::layer_thickness_variable_name, "", m_sys) - 1;

      initialize(regrid_file, (int)last_record, true);

      details::renormalize(ice_thickness, *m_layer_thickness);

      return;
    }

    auto time = m_grid->ctx()->time();

    auto requested_times = deposition_times(*m_config, *time);

    auto N_bootstrap        = static_cast<int>(m_config->get_number(N_boot_parameter));
    auto N_max              = static_cast<int>(m_config->get_number(N_max_parameter));
    auto N_deposition_times = requested_times.size();

    if (N_bootstrap < 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s have to be non-negative (got %d)",
                                    N_boot_parameter, N_bootstrap);
    }

    m_log->message(2,
                   "* Bootstrapping the isochrone tracking model, adding %d isochronal layers...\n",
                   N_bootstrap);

    if (N_bootstrap + (int)N_deposition_times > N_max) {
      auto times = m_config->get_string(times_parameter);
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "%s (%d) + %s (%d) exceeds the amount of storage allocated (%s = %d)", N_boot_parameter,
          (int)N_bootstrap, times_parameter, (int)N_deposition_times, N_max_parameter, (int)N_max);
    }

    if (N_bootstrap > 0) {
      // prepend "bootstrapping" layers to m_deposition_times
      double T_0 = time->start();

      m_top_layer_index = 0;
      std::vector<double> deposition_times;
      for (int k = 0; k < N_bootstrap; ++k) {
        deposition_times.push_back(T_0);
        m_top_layer_index += 1;
      }
      for (const auto &t : requested_times) {
        deposition_times.push_back(t);
      }

      // re-allocate storage
      m_layer_thickness = allocate_layer_thickness(m_grid, deposition_times);
      m_tmp             = m_layer_thickness->duplicate(array::WITH_GHOSTS);

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

      m_layer_thickness = allocate_layer_thickness(m_grid, requested_times);
      m_tmp             = m_layer_thickness->duplicate(array::WITH_GHOSTS);
      m_top_layer_index = 0;

      array::AccessScope scope{ &ice_thickness, m_layer_thickness.get() };

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_layer_thickness->get_column(i, j)[0] = ice_thickness(i, j);
      }
    }

  } catch (RuntimeError &e) {
    e.add_context("bootstrapping the isochrone tracking model");
    throw;
  }
}

/*!
 * Initialize from the `input_file` using `record`. Uses regridding if `use_interpolation`
 * is true.
 */
void Isochrones::initialize(const File &input_file, int record, bool use_interpolation) {
  try {
    using namespace details;

    m_log->message(2, "* Initializing the isochrone tracking model from '%s'...\n",
                   input_file.filename().c_str());

    if (use_interpolation) {
      m_log->message(2, "  [using bilinear interpolation to read layer thicknesses]\n");
    }

    const auto &time = m_grid->ctx()->time();

    {
      std::shared_ptr<array::Array3D> tmp;

      if (use_interpolation) {
        tmp = regrid_layer_thickness(m_grid, input_file, record);
      } else {
        tmp = read_layer_thickness(m_grid, input_file, record);
      }

      m_layer_thickness =
          allocate_layer_thickness(*tmp, time->start(), deposition_times(*m_config, *time));
    }
    m_tmp = m_layer_thickness->duplicate(array::WITH_GHOSTS);

    // set m_top_layer_index
    m_top_layer_index = n_active_layers(m_layer_thickness->get_levels(), time->start()) - 1;

  } catch (RuntimeError &e) {
    e.add_context("initializing the isochrone tracking model");
    throw;
  }
}

/*!
 * Re-start the model from a PISM output file.
 */
void Isochrones::restart(const File &input_file, int record) {
  if (details::regridp(*m_config)) {
    File regrid_file(m_grid->com, m_config->get_string("input.regrid.file"), io::PISM_GUESS,
                     io::PISM_READONLY);

    auto last_record = regrid_file.nrecords(details::layer_thickness_variable_name, "", m_sys);

    initialize(regrid_file, (int)last_record, true);
  } else {
    initialize(input_file, record, false);
  }
}

/*!
 * Update layer thicknesses.
 *
 * @param[in] t time step start, seconds
 * @param[in] dt time step length, seconds
 * @param[in] u x-component of the ice velocity, m/s
 * @param[in] v y-component of the ice velocity, m/s
 * @param[in] ice thickness, m
 * @param[in] top_surface_mass_balance total top surface mass balance over the time step, meters
 * @param[in] bottom_surface_mass_balance total bottom surface mass balance over the time step, meters
 *
 * For mass balance inputs, positive corresponds to mass gain.
 */
void Isochrones::update(double t, double dt, const array::Array3D &u, const array::Array3D &v,
                        const array::Scalar &ice_thickness,
                        const array::Scalar &top_surface_mass_balance,
                        const array::Scalar &bottom_surface_mass_balance) {

  // apply top surface and basal mass balance terms:
  {
    array::AccessScope scope{ &top_surface_mass_balance, &bottom_surface_mass_balance,
                              m_layer_thickness.get() };

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double *H = m_layer_thickness->get_column(i, j);

      // apply the surface mass balance
      {
        double dH = top_surface_mass_balance(i, j);

        // apply thickness change to a layer, starting from the top-most
        for (int k = (int)m_top_layer_index; k >= 0; --k) {
          if (H[k] + dH >= 0.0) {
            // thickness change is non-negative or does not remove the whole layer: apply to
            // the current layer and stop
            H[k] += dH;
            break;
          }

          dH += H[k];
          H[k] = 0.0;
        }
      }
      // apply the basal melt rate
      {
        double dH = bottom_surface_mass_balance(i, j);

        // apply thickness change to a layer, starting from the bottom
        for (size_t k = 0; k <= m_top_layer_index; ++k) {
          if (H[k] + dH >= 0.0) {
            // thickness change is non-negative or does not remove the whole layer: apply to
            // the current layer and stop
            H[k] += dH;
            break;
          }

          dH += H[k];
          H[k] = 0.0;
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
    double T                     = t + dt;
    const auto &deposition_times = m_layer_thickness->get_levels();
    size_t N                     = deposition_times.size();

    // Find the index k such that deposition_times[k] <= T < deposition_times[k + 1]
    //
    // Note: `k` below will be strictly less than `N - 1`, ensuring that the index "k + 1"
    // is valid.
    //
    // FIXME: consider using a gsl_interp_accel to speed this up
    size_t k = gsl_interp_bsearch(deposition_times.data(), T, 0, N - 1);

    double T_k = deposition_times[k];

    double top_layer_deposition_time = deposition_times.at(m_top_layer_index);
    if (T_k > top_layer_deposition_time) {
      // we reached the next requested deposition time

      if (m_top_layer_index < N - 1) {
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
  return std::min(max_timestep_deposition_times(t), max_timestep_cfl());
}

MaxTimestep Isochrones::max_timestep_cfl() const {

  if (m_stress_balance == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Isochrone tracking: no stress balance provided. "
                                  "Cannot compute the maximum time step.");
  }

  return MaxTimestep(m_stress_balance->max_timestep_cfl_3d().dt_max.value(), "isochrones");
}

/*!
 * Maximum time step we can take at time `t`.
 *
 * We can go up to the next deposition time.
 */
MaxTimestep Isochrones::max_timestep_deposition_times(double t) const {
  auto deposition_times = m_layer_thickness->get_levels();

  double t0 = deposition_times[0];
  if (t < t0) {
    return { t0 - t, "isochrones" };
  }

  if (t >= deposition_times.back()) {
    return { "isochrones" };
  }

  auto N = deposition_times.size();

  // Find the index k such that m_deposition_times[k] <= T < m_deposition_times[k + 1]
  //
  // Note: `k` below will be strictly less than `N - 1`, ensuring that the index "k + 1"
  // is valid.
  //
  // FIXME: consider using gsl_interp_accel to speed this up
  size_t k = gsl_interp_bsearch(deposition_times.data(), t, 0, N - 1);

  return { deposition_times[k + 1] - t, "isochrones" };
}

/*!
 * Define the model state in an output file.
 *
 * We are saving layer thicknesses, deposition times, and the number of active layers.
 */
void Isochrones::define_model_state_impl(const File &output) const {
  m_layer_thickness->define(output, io::PISM_DOUBLE);
}

/*!
 * Write the model state to an output file.
 */
void Isochrones::write_model_state_impl(const File &output) const {
  m_layer_thickness->write(output);
}

const array::Array3D &Isochrones::layer_thicknesses() const {
  return *m_layer_thickness;
}

namespace diagnostics {

/*! @brief Report isochrone depth below surface, in meters */
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
             Diagnostic::Ptr(new diagnostics::IsochroneDepths(this)) },
           { details::layer_thickness_variable_name, Diagnostic::wrap(*m_layer_thickness) } };
}

} // end of namespace pism

/*
 * LocalWords: LocalWords deposition
*/
