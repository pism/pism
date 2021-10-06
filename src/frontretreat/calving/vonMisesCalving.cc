/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/IceModelVec2V.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace calving {

vonMisesCalving::vonMisesCalving(IceGrid::ConstPtr grid,
                                 std::shared_ptr<const rheology::FlowLaw> flow_law)
  : StressCalving(grid, 2),
    m_calving_threshold(m_grid, "vonmises_calving_threshold", WITHOUT_GHOSTS),
    m_flow_law(flow_law)
{

  if (m_config->get_flag("calving.vonmises_calving.use_custom_flow_law")) {
    EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();
    rheology::FlowLawFactory factory("calving.vonmises_calving.", m_config, EC);
    m_flow_law = factory.create();
  }

  m_calving_rate.metadata().set_name("vonmises_calving_rate");
  m_calving_rate.set_attrs("diagnostic",
                           "horizontal calving rate due to von Mises calving",
                           "m s-1", "m year-1", "", 0);

  m_calving_threshold.set_attrs("diagnostic",
                                "threshold used by the 'von Mises' calving method",
                                "Pa", "Pa",
                                "", 0); // no standard name
  m_calving_threshold.set_time_independent(true);

}

vonMisesCalving::~vonMisesCalving() {
  // empty
}

void vonMisesCalving::init() {

  m_log->message(2,
                 "* Initializing the 'von Mises calving' mechanism...\n");

  std::string threshold_file = m_config->get_string("calving.vonmises_calving.threshold_file");
  double sigma_max = m_config->get_number("calving.vonmises_calving.sigma_max");

  m_calving_threshold.set(sigma_max);

  if (not threshold_file.empty()) {
    m_log->message(2,
                   "  Reading von Mises calving threshold from file '%s'...\n",
                   threshold_file.c_str());

    m_calving_threshold.regrid(threshold_file, CRITICAL);
  } else {
    m_log->message(2,
                   "  von Mises calving threshold: %3.3f Pa.\n", sigma_max);
  }

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "-calving vonmises_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

  m_strain_rates.set(0.0);
}

//! \brief Uses principal strain rates to apply the "von Mises" calving method
/*!
  See equation (4) in [@ref Morlighem2016].
*/
void vonMisesCalving::update(const IceModelVec2CellType &cell_type,
                             const IceModelVec2S &ice_thickness,
                             const IceModelVec2V &ice_velocity,
                             const IceModelVec3 &ice_enthalpy) {

  using std::max;
  using std::sqrt;
  using std::pow;

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  // make a copy with a wider stencil
  m_cell_type.copy_from(cell_type);

  stressbalance::compute_2D_principal_strain_rates(ice_velocity, m_cell_type,
                                                   m_strain_rates);
  m_strain_rates.update_ghosts();

  IceModelVec::AccessList list{&ice_enthalpy, &ice_thickness, &m_cell_type, &ice_velocity,
                               &m_strain_rates, &m_calving_rate, &m_calving_threshold};

  const double *z = &m_grid->z()[0];

  double glen_exponent = m_flow_law->exponent();

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (m_cell_type.ice_free_ocean(i, j) and m_cell_type.next_to_ice(i, j)) {

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
          if (m_cell_type.icy(I, j)) {
            velocity_magnitude += ice_velocity(I, j).magnitude();
            {
              double H = ice_thickness(I, j);
              unsigned int k = m_grid->kBelowHeight(H);
              hardness += averaged_hardness(*m_flow_law, H, k, &z[0], ice_enthalpy.get_column(I, j));
            }
            eigen1 += m_strain_rates(I, j, 0);
            eigen2 += m_strain_rates(I, j, 1);
            N += 1;
          }
        }

        for (int q = -1; q < 2; q += 2) {
          const int J = j + q * offset;
          if (m_cell_type.icy(i, J)) {
            velocity_magnitude += ice_velocity(i, J).magnitude();
            {
              double H = ice_thickness(i, J);
              unsigned int k = m_grid->kBelowHeight(H);
              hardness += averaged_hardness(*m_flow_law, H, k, &z[0], ice_enthalpy.get_column(i, J));
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
      const double effective_tensile_strain_rate = sqrt(0.5 * (pow(max(0.0, eigen1), 2) +
                                                               pow(max(0.0, eigen2), 2)));
      // [\ref Morlighem2016] equation 7
      const double sigma_tilde = sqrt(3.0) * hardness * pow(effective_tensile_strain_rate,
                                                            1.0 / glen_exponent);

      // Calving law [\ref Morlighem2016] equation 4
      m_calving_rate(i, j) = velocity_magnitude * sigma_tilde / m_calving_threshold(i, j);

    } else { // end of "if (ice_free_ocean and next_to_floating)"
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points
}

const IceModelVec2S& vonMisesCalving::threshold() const {
  return m_calving_threshold;
}

DiagnosticList vonMisesCalving::diagnostics_impl() const {
  return {{"vonmises_calving_rate", Diagnostic::wrap(m_calving_rate)},
          {"vonmises_calving_threshold", Diagnostic::wrap(m_calving_threshold)}};
}

} // end of namespace calving
} // end of namespace pism
