/* Copyright (C) 2016 PISM Authors
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

#include "vonMisesCalving.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceModelVec2CellType.hh"
#include "remove_narrow_tongues.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/util/PISMVars.hh"
#include "base/rheology/FlowLaw.hh"

namespace pism {
namespace calving {

vonMisesCalving::vonMisesCalving(IceGrid::ConstPtr g,
                           stressbalance::StressBalance *stress_balance)
  : StressCalving(g, stress_balance, 2) {

  m_sigma_max = m_config->get_double("calving.vonmises.sigma_max");
}

vonMisesCalving::~vonMisesCalving() {
  // empty
}

void vonMisesCalving::init() {

  m_log->message(2,
                 "* Initializing the 'von Mises calving' mechanism...\n");

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "-calving vonmises_calving using a non-square grid cell is not implemented (yet);\n"
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
void vonMisesCalving::compute_calving_rate(const IceModelVec2CellType &mask,
                                           IceModelVec2S &result) const {

  using std::max;

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  update_strain_rates();

  const IceModelVec2V &ssa_velocity  = m_stress_balance->advective_velocity();
  const IceModelVec3  *enthalpy      = m_grid->variables().get_3d_scalar("enthalpy");
  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*enthalpy);
  list.add(ice_thickness);
  list.add(mask);
  list.add(ssa_velocity);
  list.add(m_strain_rates);
  list.add(result);

  const double *z = &m_grid->z()[0];
  const rheology::FlowLaw*
    flow_law = m_stress_balance->shallow()->flow_law();

  const double ssa_n = flow_law->exponent();

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) and mask.next_to_floating_ice(i, j)) {

      double
        velocity_magnitude = 0.0,
        hardness           = 0.0;
      // Average of strain-rate eigenvalues in adjacent floating grid cells.
      double
        eigen1             = 0.0,
        eigen2             = 0.0;
      {
        int N = 0;
        for (int p = -1; p < 2; p += 2) {
          const int I = i + p * offset;
          if (mask.floating_ice(I, j)) {
            velocity_magnitude += ssa_velocity(I, j).magnitude();
            {
              double H = ice_thickness(I, j);
              unsigned int k = m_grid->kBelowHeight(H);
              hardness += averaged_hardness(*flow_law, H, k, &z[0], enthalpy->get_column(I, j));
            }
            eigen1 += m_strain_rates(I, j, 0);
            eigen2 += m_strain_rates(I, j, 1);
            N += 1;
          }
        }

        for (int q = -1; q < 2; q += 2) {
          const int J = j + q * offset;
          if (mask.floating_ice(i, J)) {
            velocity_magnitude += ssa_velocity(i, J).magnitude();
            {
              double H = ice_thickness(i, J);
              unsigned int k = m_grid->kBelowHeight(H);
              hardness += averaged_hardness(*flow_law, H, k, &z[0], enthalpy->get_column(i, J));
            }
            eigen1 += m_strain_rates(i, J, 0);
            eigen2 += m_strain_rates(i, J, 1);
            N += 1;
          }
        }

        if (N > 0) {
          eigen1             /= N;
          eigen2             /= N;
          hardness           /= N;
          velocity_magnitude /= N;
        }
      }

      // [\ref Morlighem2016] equation 6
      const double effective_tensile_strain_rate = sqrt(0.5 * (PetscSqr(max(0.0, eigen1)) +
                                                               PetscSqr(max(0.0, eigen2))));
      // [\ref Morlighem2016] equation 7
      const double sigma_tilde = sqrt(3.0) * hardness * pow(effective_tensile_strain_rate,
                                                            1.0 / ssa_n);

      // Calving law [\ref Morlighem2016] equation 4
      result(i, j) = velocity_magnitude * sigma_tilde / m_sigma_max;

    } else { // end of "if (ice_free_ocean and next_to_floating)"
      result(i, j) = 0.0;
    }
  }   // end of loop over grid points
}

std::map<std::string, Diagnostic::Ptr> vonMisesCalving::diagnostics_impl() const {
  return {{"vonmises_calving_rate",
        Diagnostic::Ptr(new CalvingRate(this, "vonmises_calving_rate",
                                        "horizontal calving rate due to von Mises calving"))}};
}

} // end of namespace calving
} // end of namespace pism
