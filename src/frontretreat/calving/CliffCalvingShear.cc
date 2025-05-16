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

#include "pism/frontretreat/calving/CliffCalvingShear.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/array/CellType.hh"

namespace pism {
namespace calving {

CliffCalvingShear::CliffCalvingShear(std::shared_ptr<const Grid> grid)
  : Component(grid),
    m_calving_rate(grid, "shear_cliff_calving_rate")  // where else do I have to change this?
{
  m_calving_rate.metadata(0)
      .long_name("horizontal calving rate due to shear stress failure")
      .units("m s^-1")
      .output_units("m year^-1");

m_C0.metadata(0)
      .long_name("scaling factor for calving rate due to shear stress failure")
      .units("m s^-1")
      .output_units("m year^-1");

m_max_cliff_calving_rate.metadata(0)
      .long_name("maximum cliff calving rate due to mélange buttressing")
      .units("m s^-1")
      .output_units("m year^-1");
}

void CliffCalvingShear::init() {

  m_log->message(2,
                 "* Initializing the 'Shear stress cliff calving' mechanism...\n");

  m_C0 = m_config->get_number("calving.cliff_calving_shear.C0");
  m_max_cliff_calving_rate = m_config->get_number("calving.cliff_calving_shear.max_cliff_calving_rate");

  m_log->message(2,
                 "  Scaling factor C0: %3.3f m/yr.\n", 
                 convert(m_sys, m_C0, "m second-1", "m year-1"));
  m_log->message(2,
                 "  Maximum cliff calving rate: %3.3f m/yr.\n",
                 convert(m_sys, m_max_cliff_calving_rate, "m second-1", "m year-1"));

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "-calving cliff_calving_shear using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

}

void CliffCalvingShear::update(const array::CellType1 &cell_type,
                             const array::Scalar &ice_thickness,
                             const array::Scalar &sea_level,
                             const array::Scalar &bed_elevation) {

  using std::min;

  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    gravity       = m_config->get_number("constants.standard_gravity");

  array::AccessScope list{&ice_thickness, &cell_type, &m_calving_rate, &sea_level,
                               &bed_elevation, &surface_elevation};

  for (auto pt = m_grid->points(); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have grounded ice neighbors after the mass continuity step
    if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_grounded_ice(i, j)) {
      // Get the ice thickness, surface elevation, and mask in all neighboring grid cells
      stencils::Star<double> H = ice_thickness.star(i, j);
      stencils::Star<double> h = surface_elevation.star(i, j);
      stencils::Star<int> M = cell_type.star_int(i, j);

      // Get the ice thickness and mask in the partially filled grid cell where we apply calving
      // it is calculated as the average of the ice thickness and surface elevation of the adjacent icy cells 
      const double H_threshold = part_grid_threshold_thickness(M, H, h, bed_elevation(i, j));
      const int m = mask::grounded_ice(sea_level(i, j), bed_elevation(i, j), H_threshold);
      
      // Calculate the parameters for the calving law given in [\ref Schlemm2019]
      const double F = H_threshold - (sea_level(i, j) - bed_elevation(i, j)),  
                  w = (sea_level(i, j) - bed_elevation(i, j)) / H_threshold,       
                  Fs = 115 * pow(w-0.356, 4) + 21,
                  Fc = 75-49*w,
                  s = 0.17 * pow(9.1, w) + 1.76,
                  Fscaled = std::max(0.0,F-Fc)/Fs;  

      // Calculate the calving rate [\ref Schlemm2019] if cell is grounded
      double C_unbuttressed;
      C_unbuttressed = (mask::grounded(m) ?
                         m_C0 * pow(Fscaled, s):
                         0.0);
                         
      // Apply mélange buttressing [\ref Schlemm2021]
      // check value of max_cliff_calving to prevent division by zero 
      m_calving_rate(i, j) =  (std::abs(m_max_cliff_calving_rate) < 1e-12 ?
                        0.0 :
                        C_unbuttressed / (1.0 + C_unbuttressed / m_max_cliff_calving_rate));
    } else {
      m_calving_rate(i, j) = 0.0;
    }
  }   // end of loop over grid points


}

const array::Scalar &CliffCalvingShear::calving_rate() const {
  return m_calving_rate;
}

DiagnosticList CliffCalvingShear::diagnostics_impl() const {
  return {{"cliff_calving_shear_rate", Diagnostic::wrap(m_calving_rate)}};
}

} // end of namespace calving
} // end of namespace pism
