/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023 PISM Authors
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

#include "pism/frontretreat/calving/EigenCalving.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/array/CellType.hh"
#include "pism/stressbalance/StressBalance.hh"

namespace pism {
namespace calving {

EigenCalving::EigenCalving(std::shared_ptr<const Grid> grid)
  : StressCalving(grid, 2) {

  m_K = m_config->get_number("calving.eigen_calving.K");

  m_calving_rate.metadata().set_name("eigen_calving_rate");
  m_calving_rate.metadata(0)
      .long_name("horizontal calving rate due to eigen-calving")
      .units("m s^-1")
      .output_units("m year^-1");
}

void EigenCalving::init() {

  m_log->message(2, "* Initializing the 'eigen-calving' mechanism...\n");

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "-calving eigen_calving using a non-square grid cell is not implemented (yet);\n"
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
void EigenCalving::update(const array::CellType &cell_type,
                          const array::Vector1 &ice_velocity) {

  // make a copy with a wider stencil
  m_cell_type.copy_from(cell_type);
  assert(m_cell_type.stencil_width() >= 2);

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  // eigenCalvOffset allows adjusting the transition from compressive to extensive flow
  // regime
  const double eigenCalvOffset = 0.0;

  stressbalance::compute_2D_principal_strain_rates(ice_velocity, m_cell_type,
                                                   m_strain_rates);
  m_strain_rates.update_ghosts();

  array::AccessScope list{&m_cell_type, &m_calving_rate, &m_strain_rates};

  // Compute the horizontal calving rate
  for (auto pt = m_grid->points(); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (m_cell_type.ice_free_ocean(i, j) and m_cell_type.next_to_floating_ice(i, j)) {

      // Average of strain-rate eigenvalues in adjacent floating grid cells to be used for
      // eigen-calving:
      double
        eigen1 = 0.0,
        eigen2 = 0.0;
      {
        int N = 0;
        for (int p = -1; p < 2; p += 2) {
          const int I = i + p * offset;
          if (m_cell_type.floating_ice(I, j) and not m_cell_type.ice_margin(I, j)) {
            eigen1 += m_strain_rates(I, j).eigen1;
            eigen2 += m_strain_rates(I, j).eigen2;
            N += 1;
          }
        }

        for (int q = -1; q < 2; q += 2) {
          const int J = j + q * offset;
          if (m_cell_type.floating_ice(i, J) and not m_cell_type.ice_margin(i, J)) {
            eigen1 += m_strain_rates(i, J).eigen1;
            eigen2 += m_strain_rates(i, J).eigen2;
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
        m_calving_rate(i, j) = m_K * eigen1 * (eigen2 - eigenCalvOffset);
      } else {
        m_calving_rate(i, j) = 0.0;
      }

    } else { // end of "if (ice_free_ocean and next_to_floating)"
      m_calving_rate(i, j) = 0.0;
    }
  } // end of the loop over grid points
}

DiagnosticList EigenCalving::diagnostics_impl() const {
  return {{"eigen_calving_rate", Diagnostic::wrap(m_calving_rate)}};
}

} // end of namespace calving
} // end of namespace pism
