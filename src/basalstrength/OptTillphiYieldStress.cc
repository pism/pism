// Copyright (C) 2004--2023 PISM Authors
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

#include "pism/basalstrength/OptTillphiYieldStress.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

/*! Optimization of till friction angle for given target surface elevation, analogous to
  Pollard et al. (2012), TC 6(5), "A simple inverse method for the distribution of basal
  sliding coefficients under ice sheets, applied to Antarctica"
*/
OptTillphiYieldStress::OptTillphiYieldStress(std::shared_ptr<const Grid> grid)
  : MohrCoulombYieldStress(grid),
    m_mask(m_grid, "diff_mask"),
    m_usurf_difference(m_grid, "usurf_difference"),
    m_usurf_target(m_grid, "usurf")
{
  // In this model tillphi is NOT time-independent.
  m_till_phi.set_time_independent(false);

  m_name = "Iterative optimization of the till friction angle for the Mohr-Coulomb yield stress model";

  m_usurf_target.metadata()
      .long_name("target surface elevation for tillphi optimization")
      .units("m");

  m_usurf_target.set_time_independent(true);

  m_usurf_difference.metadata()
      .long_name("difference between modeled and target surface elevations")
      .units("m");

  m_usurf_difference.set(0.0);

  m_mask.metadata()
      .long_name(
          "one if the till friction angle was updated by the last iteration, zero otherwise ");

  m_mask.metadata()["flag_values"] = {0.0, 1.0};
  m_mask.metadata()["flag_meanings"] = "no_update updated_during_last_iteration";

  {
    // convergence threshold
    m_dhdt_min  = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.dhdt_min", "m / s");

    // scale used to compute tillphi adjustment using the surface elevation mismatch:
    m_dphi_scale = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.dphi_scale", "degree / m");
    // upper and lower bounds of the tillphi adjustment during an iteration:
    m_dphi_max = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.dphi_max");
    m_dphi_min = -2 * m_dphi_max;
    // lower bound of tillphi:
    m_phi0_min = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.phi0_min");
    m_phi0_max = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.phi0_max");
    m_topg_min = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.topg_min");
    m_topg_max = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.topg_max");
    // upper bound of tillphi:
    m_phi_max  = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.phi_max");

    if (m_phi0_min >= m_phi_max) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "basal_yield_stress.mohr_coulomb.tillphi_opt: phi0_min >= phi_max");
    }

    if (m_topg_min >= m_topg_max) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "basal_yield_stress.mohr_coulomb.tillphi_opt: topg_min >= topg_max");
    }
  }

  {
    m_time_name = m_config->get_string("time.dimension_name") + "_tillphi_opt";
    m_t_last = time().current();
    m_update_interval = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.dt", "seconds");
    m_t_eps = m_config->get_number("time_stepping.resolution", "seconds");
  }

  m_log->message(2,
                 "  Using iterative optimization of the till friction angle.\n"
                 "  Lower bound phi0 of the till friction angle is a piecewise-linear function of bed elevation (b):\n"
                 "          /  %5.2f                                for b < %.f\n"
                 "   phi0 = |  %5.2f + (b - (%.f)) * (%.2f / %.f)   for %.f < b < %.f\n"
                 "          \\  %5.2f                               for %.f < b\n",
                 m_phi0_min, m_topg_min,
                 m_phi0_min, m_topg_min, m_phi0_max - m_phi0_min, m_topg_max - m_topg_min, m_topg_min, m_topg_max,
                 m_phi0_max, m_topg_max);
}

/*!
 * Initialize the last time tillphi was updated.
 */
void OptTillphiYieldStress::init_t_last(const File &input_file) {
  if (input_file.find_variable(m_time_name)) {
    input_file.read_variable(m_time_name, {0}, {1}, &m_t_last);
  } else {
    m_t_last = time().current();
  }
}

/*!
 * Initialize the target ice surface elevation.
 */
void OptTillphiYieldStress::init_usurf_target(const File &input_file) {
  auto filename = m_config->get_string("basal_yield_stress.mohr_coulomb.tillphi_opt.file");

  if (not filename.empty()) {
    m_usurf_target.regrid(filename, io::CRITICAL);
  } else {
    m_log->message(2, "* No file set to read target surface elevation from... using '%s'\n",
                   input_file.filename().c_str());

    m_usurf_target.regrid(input_file, io::CRITICAL);
  }

  m_usurf_target.metadata().set_name("usurf_target");
}

void OptTillphiYieldStress::restart_impl(const File &input_file, int record) {

  MohrCoulombYieldStress::restart_impl(input_file, record);

  init_t_last(input_file);

  init_usurf_target(input_file);
}

//! Initialize the pseudo-plastic till mechanical model.
//! Target surface elevation and initial iteration time are set.
void OptTillphiYieldStress::bootstrap_impl(const File &input_file,
                                           const YieldStressInputs &inputs) {

  MohrCoulombYieldStress::bootstrap_impl(input_file, inputs);

  init_t_last(input_file);

  init_usurf_target(input_file);
}

void OptTillphiYieldStress::init_impl(const YieldStressInputs &inputs) {
  (void) inputs;
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "not implemented: till friction angle optimization "
                                "cannot be initialized without an input file");
}

MaxTimestep OptTillphiYieldStress::max_timestep_impl(double t) const {
  MaxTimestep dt_max;
  {
    if (t < m_t_last) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "time %f is less than the previous time %f",
                                    t, m_t_last);
    }

    // Find the smallest time of the form m_t_last + k * m_update_interval that is greater
    // than t
    double k = ceil((t - m_t_last) / m_update_interval);

    double
      t_next = m_t_last + k * m_update_interval,
      dt = t_next - t;

    if (dt < m_t_eps) {
      dt = m_update_interval;
    }

    dt_max = MaxTimestep(dt, "tillphi_opt");
  }

  auto dt_mohr_coulomb = MohrCoulombYieldStress::max_timestep_impl(t);

  return std::min(dt_max, dt_mohr_coulomb);
}

void OptTillphiYieldStress::update_impl(const YieldStressInputs &inputs,
                                        double t, double dt) {

  double
    t_next  = m_t_last + m_update_interval,
    t_final = t + dt;

  if (t_final < m_t_last) {
    throw RuntimeError(PISM_ERROR_LOCATION, "cannot go back in time");
  }

  if (std::abs(t_next - t_final) < m_t_eps) { // reached the next update time
    update_tillphi(inputs.geometry->ice_surface_elevation,
                   inputs.geometry->bed_elevation,
                   inputs.geometry->cell_type);
    m_t_last = t_final;
  }

  MohrCoulombYieldStress::update_impl(inputs, t, dt);
}

//! Perform an iteration to adjust the till friction angle according to the difference
//! between target and modeled surface elevations.
void OptTillphiYieldStress::update_tillphi(const array::Scalar &ice_surface_elevation,
                                           const array::Scalar &bed_topography,
                                           const array::CellType &cell_type) {

  m_log->message(2, "* Updating till friction angle...\n");

  array::AccessScope list
    { &m_till_phi, &m_usurf_target, &m_usurf_difference, &m_mask, &ice_surface_elevation, &bed_topography, &cell_type };

  m_mask.set(0.0);

  double slope = (m_phi0_max - m_phi0_min) / (m_topg_max - m_topg_min);

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Compute the lower bound of the till friction angle (default value corresponds
    // to "bed_topography(i, j) > topg_max"):
    double phi0 = m_phi0_max;
    if (bed_topography(i, j) <= m_topg_min) {
      phi0 = m_phi0_min;
    } else if (bed_topography(i, j) <= m_topg_max) {
      phi0 = m_phi0_min + (bed_topography(i, j) - m_topg_min) * slope;
    }

    if (cell_type.grounded_ice(i, j)) {
      double dh_previous       = m_usurf_difference(i, j);
      m_usurf_difference(i, j) = m_usurf_target(i, j) - ice_surface_elevation(i, j);
      double dh_change         = std::abs(m_usurf_difference(i, j) - dh_previous);

      if (dh_change / m_update_interval > m_dhdt_min) {
        // Update tillphi if the rate of change of the surface elevation mismatch since
        // the last iteration is greater than the convergence threshold m_dhdt_min.

        // Mark this location as "updated by the last iteration":
        m_mask(i, j) = 1.0;

        // Compute (and clip) the tillphi adjustment
        double dphi = m_usurf_difference(i, j) * m_dphi_scale;
        dphi        = pism::clip(dphi, m_dphi_min, m_dphi_max);

        // Update (and clip) the till friction angle:
        m_till_phi(i, j) += dphi;
        m_till_phi(i, j) = pism::clip(m_till_phi(i, j), phi0, m_phi_max);
      }
    } else if (cell_type.ocean(i, j)) {
      // Floating and ice free ocean: use the bed-elevation-dependent lower bound of
      // tillphi:
      m_till_phi(i, j) = phi0;
    }
  } // end of the loop over grid points
}

void OptTillphiYieldStress::define_model_state_impl(const File &output) const {
  MohrCoulombYieldStress::define_model_state_impl(output);

  if (not output.find_variable(m_time_name)) {
    output.define_variable(m_time_name, io::PISM_DOUBLE, {});

    output.write_attribute(m_time_name, "long_name",
                           "time of the last update of the till friction angle");
    output.write_attribute(m_time_name, "calendar", time().calendar());
    output.write_attribute(m_time_name, "units", time().units_string());
  }
}

void OptTillphiYieldStress::write_model_state_impl(const File &output) const {
  MohrCoulombYieldStress::write_model_state_impl(output);

  output.write_variable(m_time_name, {0}, {1}, &m_t_last);
}

DiagnosticList OptTillphiYieldStress::diagnostics_impl() const {

  return combine({{"tillphi", Diagnostic::wrap(m_till_phi)},
                  {"usurf_difference", Diagnostic::wrap(m_usurf_difference)},
                  {"usurf_target", Diagnostic::wrap(m_usurf_target)},
                  {"diff_mask", Diagnostic::wrap(m_mask)}},
    YieldStress::diagnostics_impl());
}

} // end of namespace pism
