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

#include "pism/frontretreat/calving/CliffCalvingTensile.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/stencils.hh"
#include "pism/util/Mask.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"

namespace pism {
namespace calving {

CliffCalvingTensile::CliffCalvingTensile(std::shared_ptr<const Grid> grid)
  : Component(grid),
    m_calving_rate(grid, "tensile_cliff_calving_rate"),
    m_I(0.0),
    m_alpha(0.0)
{
  m_calving_rate.metadata(0)
      .long_name("horizontal calving rate due to tensile stress failure")
      .units("m s^-1")
      .output_units("m year^-1");
}

void CliffCalvingTensile::init() {

  m_log->message(2,
                 "* Initializing the 'Tensile stress cliff calving' mechanism...\n");

  m_I = m_config->get_number("calving.cliff_calving_tensile.I");
  m_alpha = m_config->get_number("calving.cliff_calving_tensile.alpha");

  m_log->message(2,
                 "  Scaling factor I: %3.3e m/day.\n", 
                 convert(m_sys, m_I, "m second-1", "m day-1"));
  m_log->message(2,
                 "  Exponent: %3.3f.\n",
                 m_alpha);

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "-calving cliff_calving_tensile using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }
}

void CliffCalvingTensile::update(const array::CellType1 &cell_type,
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
      // Calculate the calving rate [\ref Crawford2021] if cell is grounded
      // and cliff height is greater than 135 m
      m_calving_rate(i, j) = (mask::grounded_ice(m) && Hc > 135.0 ?
                         m_I * pow(Hc, m_alpha):
                         0.0);
                         
    } else {
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points
}

const array::Scalar &CliffCalvingTensile::calving_rate() const {
  return m_calving_rate;
}

DiagnosticList CliffCalvingTensile::diagnostics_impl() const {
  return {{"cliff_calving_tensile_rate", Diagnostic::wrap(m_calving_rate)}};
}

} // end of namespace calving
} // end of namespace pism
