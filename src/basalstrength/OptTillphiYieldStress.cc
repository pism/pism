// Copyright (C) 2004--2021 PISM Authors
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

#include "OptTillphiYieldStress.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/util/Context.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/IceModelVec2CellType.hh"
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
OptTillphiYieldStress::OptTillphiYieldStress(IceGrid::ConstPtr grid)
  : MohrCoulombYieldStress(grid),
    m_mask(m_grid, "diff_mask", WITH_GHOSTS),
    m_usurf_difference(m_grid, "usurf_difference", WITH_GHOSTS),
    m_usurf_target(m_grid, "usurf", WITH_GHOSTS)
{
  // In this model tillphi is NOT time-independent.
  m_till_phi.set_time_independent(false);

  m_name = "Iterative optimization of the till friction angle for the Mohr-Coulomb yield stress model";

  m_usurf_target.set_attrs("",
                           "target surface elevation for tillphi optimization",
                           "m", "m", "" /* no standard name */, 0);
  m_usurf_target.set_time_independent(true);

  m_usurf_difference.set_attrs("diagnostic",
                               "difference between modeled and target"
                               " surface elevations",
                               "m", "m", "", 0);
  m_usurf_difference.set(0.0);

  m_mask.set_attrs("diagnostic",
                   "one if the till friction angle was"
                   " updated by the last iteration, zero otherwise ", "", "", "", 0);
  m_mask.metadata()["flag_values"] = {0.0, 1.0};
  m_mask.metadata()["flag_meanings"] = "no_update updated_during_last_iteration";

  double start_time   = m_grid->ctx()->time()->start();
  m_last_inverse_time = start_time;

  {
    // time interval between iterations:
    m_update_interval = m_config->get_number("basal_yield_stress.mohr_coulomb.tillphi_opt.dt",
                                             "seconds");

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
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -tillphi_opt arguments: phi0_min < phi_max is required");
    }

    if (m_topg_min >= m_topg_max) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -tillphi_opt arguments: topg_min < topg_max is required");
    }
  }

  m_log->message(2,
                 "  Using iterative optimization of the till friction angle.\n"
                 "  Lower bound phi0 of the till friction angle is a piecewise-linear function of bed elevation (b):\n"
                 "          /  %5.2f                                for b < %.f\n"
                 "   phi0 = |  %5.2f + (b - (%.f)) * (%.2f / %.f)   for %.f < b < %.f\n"
                 "          \\  %5.2f                               for %.f < b\n",
                 m_phi0_min, m_topg_min, m_phi0_min, m_topg_min, m_phi0_max - m_phi0_min, m_topg_max - m_topg_min, m_topg_min, m_topg_max,
                 m_phi0_max, m_topg_max);
}

void OptTillphiYieldStress::restart_impl(const File &input_file, int record) {

  MohrCoulombYieldStress::restart_impl(input_file, record);

  auto phi_file = m_config->get_string("basal_yield_stress.mohr_coulomb.tillphi_opt.file");

  if (not phi_file.empty()) {
     m_usurf_target.regrid(phi_file, CRITICAL);
  } else {
    m_log->message(2, "* No file set to read target surface elevation from... using '%s'\n",
                   input_file.filename().c_str());

    m_usurf_target.regrid(input_file, CRITICAL);
  }

  m_usurf_target.metadata().set_name("usurf_target");
}

//! Initialize the pseudo-plastic till mechanical model.
//! Target surface elevation and initial iteration time are set.
void OptTillphiYieldStress::bootstrap_impl(const File &input_file,
                                           const YieldStressInputs &inputs) {

  MohrCoulombYieldStress::bootstrap_impl(input_file, inputs);

  auto phi_file = m_config->get_string("basal_yield_stress.mohr_coulomb.tillphi_opt.file");

  if (not phi_file.empty()) {
    m_usurf_target.regrid(phi_file, CRITICAL);
  } else {
    m_log->message(2, "* No file set to read target surface elevation from... using '%s'\n",
                   input_file.filename().c_str());

    m_usurf_target.regrid(input_file, CRITICAL);
  }

  m_usurf_target.metadata().set_name("usurf_target");
}

MaxTimestep OptTillphiYieldStress::max_timestep_impl(double t) const {

  auto dt = MohrCoulombYieldStress::max_timestep_impl(t);

  if (dt.finite()) {
    return MaxTimestep(dt.value(), name());
  }

  return MaxTimestep(name());

  //! FIXME: Add max timestep according to dt_phi_inv
}

void OptTillphiYieldStress::update_impl(const YieldStressInputs &inputs,
                                        double t, double dt) {
  double dt_inverse = t - m_last_inverse_time;

  // FIXME: use predictable time step lengths
  if (dt_inverse > m_update_interval) {

    update_tillphi(inputs.geometry->ice_surface_elevation,
                       inputs.geometry->bed_elevation,
                       inputs.geometry->cell_type);

    m_last_inverse_time = t;
  }

  MohrCoulombYieldStress::update_impl(inputs, t, dt);
}

//! Perform an iteration to adjust the till friction angle according to the difference
//! between target and modeled surface elevations.
void OptTillphiYieldStress::update_tillphi(const IceModelVec2S &ice_surface_elevation,
                                           const IceModelVec2S &bed_topography,
                                           const IceModelVec2CellType &cell_type) {

  m_log->message(2, "* Updating till friction angle...\n");

  IceModelVec::AccessList list
    { &m_till_phi, &m_usurf_target, &m_usurf_difference, &m_mask, &ice_surface_elevation, &bed_topography, &cell_type };

  m_mask.set(0.0);

  double slope = (m_phi0_max - m_phi0_min) / (m_topg_max - m_topg_min);

  for (Points p(*m_grid); p; p.next()) {
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

DiagnosticList OptTillphiYieldStress::diagnostics_impl() const {

  return combine({{"tillphi", Diagnostic::wrap(m_till_phi)},
                  {"usurf_difference", Diagnostic::wrap(m_usurf_difference)},
                  {"usurf_target", Diagnostic::wrap(m_usurf_target)},
                  {"diff_mask", Diagnostic::wrap(m_mask)}},
    YieldStress::diagnostics_impl());
}

} // end of namespace pism
