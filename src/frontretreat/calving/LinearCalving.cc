/* Copyright (C) 2016, 2017, 2018, 2019, 2021, 2022, 2023 PISM Authors
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
#include "pism/util/Mask.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"

namespace pism {
namespace calving {

LinearCalving::LinearCalving(std::shared_ptr<const Grid> grid)
  : Component(grid),
    m_calving_rate(grid, "linear_calving_rate"),
    m_a(0.0),
    m_b(0.0)

{
  m_calving_rate.metadata(0)
      .long_name("horizontal calving rate with a linear dependence on the cliff height")
      .units("m s^-1")
      .output_units("m year^-1");
}

void LinearCalving::init() {

  m_log->message(2,
                 "* Initializing the 'Linear calving' mechanism...\n");

  m_a = m_config->get_number("calving.linear_calving.a");
  m_b = m_config->get_number("calving.linear_calving.b");

  // Convert from m/year to m/s
  m_b = convert(m_sys, m_b, "m year-1", "m second-1");
  m_a = convert(m_sys, m_a, "year-1", "second-1");

  m_log->message(2,
                 "  Coefficient a: %3.3f year^-1.\n", 
                 convert(m_sys, m_a, "second-1", "year-1"));
  m_log->message(2,
                 "  Constant b: %3.3f m year^-1.\n",
                 convert(m_sys, m_b, "m second-1", "m year-1"));

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "-calving linear_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
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
void LinearCalving::update(const array::CellType1 &cell_type,
                             const array::Scalar &ice_thickness,
                             const array::Scalar &sea_level,
                             const array::Scalar &bed_elevation) {

  using std::min;

  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    gravity       = m_config->get_number("constants.standard_gravity");

  GeometryCalculator gc(*m_config);

  array::AccessScope list{&ice_thickness, &cell_type, &m_calving_rate, &sea_level,
                               &bed_elevation};

  // Add statistics tracking
  double max_calving_rate = 0.0;
  double max_cliff_height = 0.0;
  int num_calving_cells = 0;
  int num_extreme_calving = 0;
  int max_rate_i = 0, max_rate_j = 0;

  for (auto pt = m_grid->points(); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have grounded ice neighbors after the mass continuity step
    if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_grounded_ice(i, j)) {
      // Get the ice thickness, surface elevation, and mask in all neighboring grid cells
      stencils::Star<double> H;   // TO DO: change H to type Scalar1 (array::Scalar1) then we can use the star stencil directly, but need  to change this also in the call to the calving update function!!!
      H.c = ice_thickness(i, j);
      H.e = ice_thickness(i+1, j);
      H.w = ice_thickness(i-1, j);
      H.n = ice_thickness(i, j+1);
      H.s = ice_thickness(i, j-1);
      stencils::Star<double> surface_elevation;
      for (auto d : {North, East, South, West}) {
        surface_elevation[d] = H[d] + bed_elevation(i, j);
      }
      stencils::Star<int> M = cell_type.star_int(i, j);

      // Get the ice thickness and mask in the partially filled grid cell where we apply calving
      // it is calculated as the average of the ice thickness and surface elevation of the adjacent icy cells 
      const double H_threshold = part_grid_threshold_thickness(M, H, surface_elevation, bed_elevation(i, j));
      const int m = gc.mask(sea_level(i, j), bed_elevation(i, j), H_threshold);
      //Cliff height
      const double Hc = H_threshold - (sea_level(i, j) - bed_elevation(i, j));
      // Calculate the calving rate [\ref Parsons2025] if cell is grounded
      m_calving_rate(i, j) = (mask::grounded_ice(m)  ?
                         std::max(0.0, m_a * Hc + m_b):
                         0.0);

      // Track statistics                   
      if (m_calving_rate(i, j) > 0.0) {
        num_calving_cells++;
        if (m_calving_rate(i, j) > max_calving_rate) {
          max_calving_rate = m_calving_rate(i, j);
          max_cliff_height = Hc;
          max_rate_i = i;
          max_rate_j = j;
        }
        // Log very high calving rates that might cause instability
        if (m_calving_rate(i, j) > 1e-5) {  // More than ~315 m/year
          num_extreme_calving++;
          m_log->message(3,
                     "! High linear calving rate at (i,j) = (%d,%d): %.2f m/year (H=%.1f m)\n",
                     i, j, m_calving_rate(i, j) * 31557600.0, Hc);
        }
      }
                         
    } else {
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points

  // Print summary statistics
  if (num_calving_cells > 0) {
    m_log->message(3,
                 "* Linear calving summary:\n"
                 "  - Active calving cells: %d\n"
                 "  - Cells with extreme rates (>315 m/year): %d\n"
                 "  - Maximum rate: %.2f m/year at (i,j)=(%d,%d)\n"
                 "  - Maximum cliff height: %.1f m\n",
                 num_calving_cells,
                 num_extreme_calving,
                 max_calving_rate * 31557600.0,
                 max_rate_i, max_rate_j,
                 max_cliff_height);
  } else {
    m_log->message(3, "* No active linear calving cells at this time step (maximum cliff height: %.1f m).\n",
                   max_cliff_height);
  }
}

const array::Scalar &LinearCalving::calving_rate() const {
  return m_calving_rate;
}

DiagnosticList LinearCalving::diagnostics_impl() const {
  return {{"linear_calving_rate", Diagnostic::wrap(m_calving_rate)}};
}

} // end of namespace calving
} // end of namespace pism
