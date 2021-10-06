/* Copyright (C) 2016, 2017, 2018, 2019, 2021 PISM Authors
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

#include "HayhurstCalving.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace calving {

HayhurstCalving::HayhurstCalving(IceGrid::ConstPtr grid)
  : Component(grid),
    m_calving_rate(grid, "hayhurst_calving_rate", WITH_GHOSTS)
{
  m_calving_rate.set_attrs("diagnostic",
                           "horizontal calving rate due to Hayhurst calving",
                           "m s-1", "m day-1", "", 0);

}

HayhurstCalving::~HayhurstCalving() {
  // empty
}

void HayhurstCalving::init() {

  m_log->message(2,
                 "* Initializing the 'Hayhurst calving' mechanism...\n");

  m_B_tilde = m_config->get_number("calving.hayhurst_calving.B_tilde");
  m_exponent_r = m_config->get_number("calving.hayhurst_calving.exponent_r");
  m_sigma_threshold = m_config->get_number("calving.hayhurst_calving.sigma_threshold", "Pa");

  m_log->message(2,
                 "  B tilde parameter: %3.3f MPa-%3.3f yr-1.\n", m_B_tilde, m_exponent_r);
  m_log->message(2,
                 "  Hayhurst calving threshold: %3.3f MPa.\n",
                 convert(m_sys, m_sigma_threshold, "Pa", "MPa"));

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "-calving hayhurst_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

}

void HayhurstCalving::update(const IceModelVec2CellType &cell_type,
                             const IceModelVec2S &ice_thickness,
                             const IceModelVec2S &sea_level,
                             const IceModelVec2S &bed_elevation) {

  using std::min;

  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    gravity       = m_config->get_number("constants.standard_gravity"),
    // convert "Pa" to "MPa" and "m yr-1" to "m s-1"
    unit_scaling  = pow(1e-6, m_exponent_r) * convert(m_sys, 1.0, "m year-1", "m second-1");

  IceModelVec::AccessList list{&ice_thickness, &cell_type, &m_calving_rate, &sea_level,
                               &bed_elevation};

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    double water_depth = sea_level(i, j) - bed_elevation(i, j);

    if (cell_type.icy(i, j) and water_depth > 0.0) {
      // note that ice_thickness > 0 at icy locations
      assert(ice_thickness(i, j) > 0);

      double H = ice_thickness(i, j);

      // Note that for ice at floatation water_depth = H * (ice_density / water_density),
      // so omega cannot exceed ice_density / water_density.
      double omega = water_depth / H;

      // Extend the calving parameterization to ice shelves. This tweak should (I hope)
      // prevent a feedback in which the creation of an ice shelf reduces the calving
      // rate, which leads to an advance of the front and an even lower calving rate, etc.
      if (omega > ice_density / water_density) {
        // ice at the front is floating
        double freeboard = (1.0 - ice_density / water_density) * H;
        H = water_depth + freeboard;
        omega = water_depth / H;
      }

      // [\ref Mercenier2018] maximum tensile stress approximation
      double sigma_0 = (0.4 - 0.45 * pow(omega - 0.065, 2.0)) * ice_density * gravity * H;

      // ensure that sigma_0 - m_sigma_threshold >= 0
      sigma_0 = std::max(sigma_0, m_sigma_threshold);

      // [\ref Mercenier2018] equation 22
      m_calving_rate(i, j) = (m_B_tilde * unit_scaling *
                              (1.0 - pow(omega, 2.8)) *
                              pow(sigma_0 - m_sigma_threshold, m_exponent_r) * H);
    } else { // end of "if (ice_free_ocean and next_to_floating)"
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points

  // Set calving rate *near* grounded termini to the average of grounded icy
  // neighbors: front retreat code uses values at these locations (the rest is for
  // visualization).

  m_calving_rate.update_ghosts();

  const Direction dirs[] = {North, East, South, West};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free(i, j) and cell_type.next_to_ice(i, j) ) {

      auto R = m_calving_rate.star(i, j);
      auto M = cell_type.star(i, j);

      int N = 0;
      double R_sum = 0.0;
      for (int n = 0; n < 4; ++n) {
        Direction direction = dirs[n];
        if (mask::icy(M[direction])) {
          R_sum += R[direction];
          N++;
        }
      }

      if (N > 0) {
        m_calving_rate(i, j) = R_sum / N;
      }
    }
  }
}

const IceModelVec2S &HayhurstCalving::calving_rate() const {
  return m_calving_rate;
}

DiagnosticList HayhurstCalving::diagnostics_impl() const {
  return {{"hayhurst_calving_rate", Diagnostic::wrap(m_calving_rate)}};
}

} // end of namespace calving
} // end of namespace pism
