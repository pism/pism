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

  m_usurf_target.set_attrs("internal",
                           "target surface elevation for tillphi optimization",
                           "m", "m", "" /* no standard name */, 0);
  m_usurf_target.set_time_independent(true);

  m_usurf_difference.set_attrs("internal",
                               "difference between modeled and target surface elevations",
                               "m", "m", "", 0);

  m_mask.set_attrs("internal", "mask for till phi iteration", "", "", "", 0);

  double start_time   = m_grid->ctx()->time()->start();
  m_last_inverse_time = start_time;

  {
    m_dt_phi_inv = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.dt",
                                        "seconds");
    m_h_inv     = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.h_inv");
    // FIXME: check units of dhdt_conv
    m_dhdt_conv = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.dh_conv", "m / s");
    m_dphi_max  = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.dphi");
    m_dphi_min  = -0.5 * m_dphi_max;
    m_phi_min   = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_min");
    m_phi_minup = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_minup");
    m_phi_max   = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_max");
    m_topg_min  = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.topg_min");
    m_topg_max  = m_config->get_number("basal_yield_stress.mohr_coulomb.iterative_phi.topg_max");

    m_slope = (m_phi_minup - m_phi_min) / (m_topg_max - m_topg_min);

    m_log->message(2,
                   "  lower bound of till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
                   "             /  %5.2f                                         for   topg < %.f\n"
                   "   phi_min = |  %5.2f + (topg - (%.f)) * (%.2f / %.f)   for   %.f < topg < %.f\n"
                   "             \\  %5.2f                                        for   %.f < topg\n",
                   m_phi_min, m_topg_min, m_phi_min, m_topg_min, m_phi_minup - m_phi_min, m_topg_max - m_topg_min, m_topg_min, m_topg_max,
                   m_phi_minup, m_topg_max);

    if (m_phi_min >= m_phi_max) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -iterative_phi arguments: phi_min < phi_max is required");
    }

    if (m_topg_min >= m_topg_max) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -iterative_phi arguments: topg_min < topg_max is required");
    }
  }
}

void OptTillphiYieldStress::restart_impl(const File &input_file, int record) {

  MohrCoulombYieldStress::restart_impl(input_file, record);

  auto phi_file = m_config->get_string("basal_yield_stress.mohr_coulomb.iterative_phi.file");

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

  auto phi_file = m_config->get_string("basal_yield_stress.mohr_coulomb.iterative_phi.file");

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
  if (dt_inverse > m_dt_phi_inv) {

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

  m_log->message(2, "\n* Perform iterative step for optimization of till friction angle phi!\n\n");

  IceModelVec::AccessList list
    { &m_till_phi, &m_usurf_target, &m_usurf_difference, &m_mask, &ice_surface_elevation, &bed_topography, &cell_type };

  m_mask.set(1.0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double usurf_difference_prev = m_usurf_difference(i, j);
    m_usurf_difference(i, j)     = ice_surface_elevation(i, j) - m_usurf_target(i, j);
    double dh_step               = std::abs(m_usurf_difference(i, j) - usurf_difference_prev);

    if (cell_type.grounded_ice(i, j)) {

      // Convergence criterion
      if (dh_step / m_dt_phi_inv > m_dhdt_conv) {
        m_mask(i, j) = 1.0;

        double dphi = m_usurf_difference(i, j) / m_h_inv;
        dphi        = pism::clip(dphi, m_dphi_min, m_dphi_max);

        m_till_phi(i, j) -= dphi;

        // default value corresponds to "bed_topography(i, j) > topg_max"
        double tillphi_min = m_phi_minup;
        if (bed_topography(i, j) <= m_topg_min) {
          tillphi_min = m_phi_min;
        } else if (bed_topography(i, j) <= m_topg_max) {
          tillphi_min = m_phi_min + (bed_topography(i, j) - m_topg_min) * m_slope;
        }

        m_till_phi(i, j) = pism::clip(m_till_phi(i, j), tillphi_min, m_phi_max);
      } else {
        m_mask(i, j) = 0.0;
      }

      // Floating and ice free ocean
    } else if (cell_type.ocean(i, j)) {
      m_till_phi(i, j) = m_phi_min;
      m_mask(i, j)     = 0.0;
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
