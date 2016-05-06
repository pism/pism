/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#include "EigenCalving.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"
#include "remove_narrow_tongues.hh"

namespace pism {
namespace calving {

EigenCalving::EigenCalving(IceGrid::ConstPtr g,
                           stressbalance::StressBalance *stress_balance)
  : CalvingFrontRetreat(g), m_stencil_width(2),
    m_stress_balance(stress_balance) {
  m_strain_rates.create(m_grid, "strain_rates", WITH_GHOSTS,
                        m_stencil_width,
                        2);

  m_strain_rates.metadata(0).set_name("eigen1");
  m_strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "second-1", "", 0);

  m_strain_rates.metadata(1).set_name("eigen2");
  m_strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "second-1", "", 1);

  m_K = m_config->get_double("eigen_calving_K");
}

EigenCalving::~EigenCalving() {
  // empty
}

void EigenCalving::init() {

  m_log->message(2,
             "* Initializing the 'eigen-calving' mechanism...\n");

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted("-calving eigen_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

  m_strain_rates.set(0.0);

}

//! \brief Uses principal strain rates to apply "eigencalving" with constant K.
/*!
  See equation (26) in [\ref Winkelmannetal2011].
*/
void EigenCalving::compute_calving_rate(const IceModelVec2CellType &mask,
                                        IceModelVec2S &result) {

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  // eigenCalvOffset allows adjusting the transition from
  // compressive to extensive flow regime
  const double eigenCalvOffset = 0.0;

  update_strain_rates();

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(result);
  list.add(m_strain_rates);

  // Compute the horizontal calving rate
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) and mask.next_to_floating_ice(i, j)) {

      // Average of strain-rate eigenvalues in adjacent floating grid cells to be used for
      // eigen-calving:
      double
        eigen1 = 0.0,
        eigen2 = 0.0;
      {
        int N = 0;
        for (int p = -1; p < 2; p += 2) {
          const int I = i + p * offset;
          if (mask.floating_ice(I, j) and not mask.ice_margin(I, j)) {
            eigen1 += m_strain_rates(I, j, 0);
            eigen2 += m_strain_rates(I, j, 1);
            N += 1;
          }
        }

        for (int q = -1; q < 2; q += 2) {
          const int J = j + q * offset;
          if (mask.floating_ice(i, J) and not mask.ice_margin(i, J)) {
            eigen1 += m_strain_rates(i, J, 0);
            eigen2 += m_strain_rates(i, J, 1);
            N += 1;
          }
        }

        if (N > 0) {
          eigen1 /= N;
          eigen2 /= N;
        }
      }

      // Calving law
      //
      // eigen1 * eigen2 has units [s^-2] and calving_rate_horizontal
      // [m*s^1] hence, eigen_calving_K has units [m*s]
      if (eigen2 > eigenCalvOffset and eigen1 > 0.0) {
        // spreading in all directions
        result(i, j) = m_K * eigen1 * (eigen2 - eigenCalvOffset);
      } else {
        result(i, j) = 0.0;
      }

    } else { // end of "if (ice_free_ocean and next_to_floating)"
      result(i, j) = 0.0;
    }
  } // end of the loop over grid points

}

void EigenCalving::add_vars_to_output_impl(const std::string &/*keyword*/, std::set<std::string> &/*result*/) {
  // empty
}

void EigenCalving::define_variables_impl(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                                  IO_Type /*nctype*/) {
  // empty
}

void EigenCalving::write_variables_impl(const std::set<std::string> &/*vars*/, const PIO& /*nc*/) {
  // empty
}

/**
 * Update the strain rates field.
 *
 * Note: this code uses the mask obtained from the Vars
 * dictionary, because the velocity field used to compute it need not
 * extend past the ice margin corresponding to the *beginning* of the
 * time-step.
 *
 * @return 0 on success
 */
void EigenCalving::update_strain_rates() {
  const IceModelVec2V        &ssa_velocity = m_stress_balance->advective_velocity();
  const IceModelVec2CellType &mask         = *m_grid->variables().get_2d_cell_type("mask");
  m_stress_balance->compute_2D_principal_strain_rates(ssa_velocity, mask, m_strain_rates);
}

} // end of namespace calving
} // end of namespace pism
