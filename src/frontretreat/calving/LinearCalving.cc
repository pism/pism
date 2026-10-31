/* Copyright (C) 2016, 2017, 2018, 2019, 2021, 2022, 2023, 2026 PISM Authors
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

#include "pism/frontretreat/calving/LinearCalving.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/stencils.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Mask.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace calving {

static const double high_calving_rate = 1e-5; // m/s

LinearCalving::LinearCalving(std::shared_ptr<const Grid> grid)
    : Component(grid), m_calving_rate(grid, "linear_calving_rate"), m_a(0.0), m_b(0.0) {

  m_use_floatation_thickness = m_config->get_flag("calving.grounded_calving.use_floatation_thickness");

  m_a = m_config->get_number("calving.linear_calving.a", "m second-1");
  m_b = m_config->get_number("calving.linear_calving.b", "m second-1");

  m_calving_rate.metadata(0)
      .long_name("horizontal calving rate with a linear dependence on the cliff height")
      .units("m s^-1")
      .output_units("m year^-1");

  double dx = m_grid->dx(), dy = m_grid->dy();
  if (fabs(dx - dy) / std::min(dx, dy) > 1e-2) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION,
        "-calving linear_calving using a non-square grid cell is not implemented (yet);\n"
        "dx = %f, dy = %f, relative difference = %f",
        dx, dy, fabs(dx - dy) / std::max(dx, dy));
  }
}

void LinearCalving::init() {

  m_log->message(2,
                 "* Initializing the 'Linear calving' mechanism...\n");
  m_log->message(2,
                 "  Coefficient a: %3.3f year^-1.\n", 
                 convert(m_sys, m_a, "second-1", "year-1"));
  m_log->message(2,
                 "  Constant b: %3.3f m year^-1.\n",
                 convert(m_sys, m_b, "m second-1", "m year-1"));

  if (m_use_floatation_thickness) {
    m_log->message(2,
                   "  Using floatation thickness for cliff height calculation.\n");
  }
}

/*!
 * @brief Update the linear calving rate
 *
 * This method computes the calving rate for each grid cell based on the
 * linear calving parameterization. The calving rate is calculated
 * only for grounded ice cells that are adjacent to ice-free ocean cells.
 *
 * The algorithm:
 * 1. Identifies ice-free ocean cells that are adjacent to grounded ice
 * 2. Calculates the cliff height above water level
 * 3. Applies the linear calving rate formula: c = a * H_c + b
 * 4. Sets negative calving rates to zero
 * 5. Tracks statistics for monitoring and debugging
 *
 * @param cell_type Ice/ocean/grounded mask
 * @param ice_thickness Ice thickness field
 * @param sea_level Sea level field
 * @param bed_elevation Bed elevation field
 */
void LinearCalving::update(const array::CellType1 &cell_type, const array::Scalar1 &ice_thickness,
                           const array::Scalar &sea_level, const array::Scalar &bed_elevation) {
  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density");

  GeometryCalculator gc(*m_config);

  array::AccessScope list{ &ice_thickness, &cell_type, &m_calving_rate, &sea_level,
                           &bed_elevation };

  // Add statistics tracking
  double max_calving_rate = 0.0;
  double max_cliff_height = 0.0;
  int num_calving_cells = 0;
  int num_extreme_calving = 0;

  for (auto pt : m_grid->points()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the ice-free ocean that have grounded
    // ice neighbors
    if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_grounded_ice(i, j)) {
      // Get the ice thickness, surface elevation, and mask in all neighboring grid cells
      auto H = ice_thickness.star(i, j);

      // at an ice-free ocean location surface elevation is equal to sea level
      stencils::Star<double> surface_elevation(sea_level(i, j));
      for (auto d : {North, East, South, West}) {
        surface_elevation[d] = H[d] + bed_elevation(i, j);
      }
      auto M = cell_type.star_int(i, j);

      // Calculate cliff height using either floatation thickness calculated from water
      // depth or threshold thickness calculated from ice thickness and surface elevation
      // of neighbouring icy cells.
      const double water_depth = sea_level(i, j) - bed_elevation(i, j);
      double cliff_height;
      int mask;
      if (m_use_floatation_thickness) {
        // Calculate floatation thickness: thickness of ice that would float in this water depth
        const double H_flotation = water_depth * (water_density / ice_density);

        // Calculate mask of the partially filled grid cell where we apply calving
        mask = gc.mask(sea_level(i, j), bed_elevation(i, j), H_flotation);

        // Calculate cliff height using floatation thickness
        cliff_height = H_flotation - water_depth;
      } else {
        // Get the ice thickness in the partially filled grid cell where we apply calving.
        // It is calculated as the average of the ice thickness and surface elevation of
        // the adjacent icy cells.
        const double H_threshold =
            part_grid_threshold_thickness(M, H, surface_elevation, bed_elevation(i, j));

        // Calculate mask of the partially filled grid cell where we apply calving
        mask = gc.mask(sea_level(i, j), bed_elevation(i, j), H_threshold);

        // Calculate cliff height (ice thickness above sea level)
        cliff_height = H_threshold - water_depth;
      }
      // Calculate the calving rate [\ref Parsons2025] if cell is grounded
      m_calving_rate(i, j) =
          (mask::grounded_ice(mask) ? std::max(0.0, m_a * cliff_height + m_b) : 0.0);

      max_cliff_height = std::max(max_cliff_height, cliff_height);
      max_calving_rate = std::max(max_calving_rate, m_calving_rate(i, j));

      // Track statistics
      num_calving_cells += m_calving_rate(i, j) > 0 ? 1 : 0;
      // Log very high (more than ~315 m/year) calving rates that might cause instability
      num_extreme_calving += (m_calving_rate(i, j) > high_calving_rate) ? 1 : 0;

    } else {
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points

  // Print summary statistics
  {
    num_calving_cells = GlobalSum(m_grid->com, num_calving_cells);
    if (num_calving_cells > 0) {
      num_extreme_calving = GlobalSum(m_grid->com, num_extreme_calving);
      max_calving_rate    = GlobalMax(m_grid->com, max_calving_rate);
      max_cliff_height    = GlobalMax(m_grid->com, max_cliff_height);
      m_log->message(3,
                     "* Linear calving summary:\n"
                     "  - Active calving cells: %d\n"
                     "  - Cells with extreme rates (>315 m/year): %d\n"
                     "  - Maximum rate: %.2f m/year\n"
                     "  - Maximum cliff height: %.1f m\n",
                     num_calving_cells, num_extreme_calving,
                     units::convert(m_sys, max_calving_rate, "m second-1", "m year-1"),
                     max_cliff_height);
    } else {
      m_log->message(
          3, "* No active linear calving cells at this time step (maximum cliff height: %.1f m).\n",
          max_cliff_height);
    }
  }
}

const array::Scalar &LinearCalving::calving_rate() const {
  return m_calving_rate;
}

DiagnosticList LinearCalving::diagnostics_impl() const {
  return { { "linear_calving_rate", Diagnostic::wrap(m_calving_rate) } };
}

} // end of namespace calving
} // end of namespace pism
