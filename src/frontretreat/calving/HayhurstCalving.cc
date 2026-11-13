/* Copyright (C) 2016, 2017, 2018, 2019, 2021, 2022, 2023, 2025, 2026 PISM Authors
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

#include "pism/frontretreat/calving/HayhurstCalving.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/Logger.hh"
#include "pism/util/stencils.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace calving {

HayhurstCalving::HayhurstCalving(std::shared_ptr<const Grid> grid)
  : Component(grid),
    m_calving_rate(grid, "hayhurst_calving_rate")
{
  m_calving_rate.metadata(0)
      .long_name("horizontal calving rate due to Hayhurst calving")
      .units("m s^-1")
      .output_units("m day^-1");
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

  
  // Read floatation thickness option
  m_use_floatation_thickness = m_config->get_flag("calving.grounded_calving.use_floatation_thickness");
  if (m_use_floatation_thickness) {
    m_log->message(2,
                   "  Using floatation thickness for cliff height calculation.\n");
  }

  double dx = m_grid->dx(), dy = m_grid->dy();
  if (fabs(dx - dy) / std::min(dx, dy) > 1e-2) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION,
        "-calving hayhurst_calving using a non-square grid cell is not implemented (yet);\n"
        "dx = %f, dy = %f, relative difference = %f",
        dx, dy, fabs(dx - dy) / std::max(dx, dy));
  }
}

void HayhurstCalving::update(const array::CellType1 &cell_type, const array::Scalar1 &ice_thickness,
                             const array::Scalar &sea_level, const array::Scalar &bed_elevation) {

  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    gravity       = m_config->get_number("constants.standard_gravity"),
    // convert "Pa" to "MPa" and "m yr-1" to "m s-1"
    unit_scaling  = pow(1e-6, m_exponent_r) * convert(m_sys, 1.0, "m year-1", "m second-1");

  array::AccessScope list{&ice_thickness, &cell_type, &m_calving_rate, &sea_level,
                          &bed_elevation};

  for (auto pt : m_grid->points()) {
    const int i = pt.i(), j = pt.j();

    double water_depth = sea_level(i, j) - bed_elevation(i, j);

    if (cell_type.icy(i, j) and water_depth > 0.0) {
      // note that ice_thickness > 0 at icy locations
      assert(ice_thickness(i, j) > 0);

      // Determine ice thickness using either floatation thickness or modelled ice thickness
      double H;
      if (m_use_floatation_thickness) {
        // Calculate floatation thickness
        const double H_flotation = water_depth * (water_density / ice_density);
        H = H_flotation;
      } else {
        H = ice_thickness(i, j);
      }


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
      // NOLINTBEGIN(readability-magic-numbers)
      double sigma_0 = (0.4 - 0.45 * pow(omega - 0.065, 2.0)) * ice_density * gravity * H;

      // ensure that sigma_0 - m_sigma_threshold >= 0
      sigma_0 = std::max(sigma_0, m_sigma_threshold);

      // [\ref Mercenier2018] equation 22
      m_calving_rate(i, j) = (m_B_tilde * unit_scaling *
                              (1.0 - pow(omega, 2.8)) *
                              pow(sigma_0 - m_sigma_threshold, m_exponent_r) * H);
      // NOLINTEND(readability-magic-numbers)

    } else { // end of "if (icy and water_depth > 0)"
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points

  // Set calving rate *near* grounded termini to the average of grounded icy
  // neighbors: front retreat code uses values at these locations (the rest is for
  // visualization).

  m_calving_rate.update_ghosts();

  // Add statistics tracking for ice-ocean boundary cells
  double max_calving_rate = 0.0;
  double max_cliff_height = 0.0;

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free(i, j) and cell_type.next_to_ice(i, j) ) {

      auto R = m_calving_rate.star(i, j);
      auto M = cell_type.star_int(i, j);
      auto H = ice_thickness.star(i, j);

      int N = 0;
      double R_sum = 0.0;
      double H_sum = 0.0;
      for (auto d : {North, East, South, West}) {
        if (mask::icy(M[d])) {
          R_sum += R[d];

          // Get ice thickness from the neighboring icy cell
          if ((d == North and j + 1 < (int)m_grid->My()) or
              (d == East  and i + 1 < (int)m_grid->Mx()) or
              (d == South and j - 1 >= 0) or
              (d == West  and i - 1 >= 0)) {
            H_sum += H[d];
          }
          N++;
        }
      }

      if (N > 0) {
        m_calving_rate(i, j) = R_sum / N;
        double avg_H = H_sum / N;

        // Track statistics for this boundary cell
        if (m_calving_rate(i, j) > max_calving_rate) {
          max_calving_rate = m_calving_rate(i, j);
          max_cliff_height = avg_H;
        }
      }
    }
  } // end of the loop over grid points

  // Print summary statistics
  {
    max_calving_rate = GlobalMax(m_grid->com, max_calving_rate);
    if (max_calving_rate > 0.0) {
      max_cliff_height = GlobalMax(m_grid->com, max_cliff_height);
      m_log->message(
          3, "* Hayhurst calving summary: max rate = %.2f m/year, max cliff height = %.1f m\n",
          units::convert(m_sys, max_calving_rate, "m second-1", "m year-1"), max_cliff_height);
    } else {
      m_log->message(3, "* No active Hayhurst calving cells at this time step.\n");
    }
  }
}

const array::Scalar &HayhurstCalving::calving_rate() const {
  return m_calving_rate;
}

DiagnosticList HayhurstCalving::spatial_diagnostics_impl() const {
  return { { "hayhurst_calving_rate", Diagnostic::wrap(m_calving_rate) } };
}

} // end of namespace calving
} // end of namespace pism
