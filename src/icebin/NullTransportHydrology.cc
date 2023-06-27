// Copyright (C) 2012-2014, 2016, 2023 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"
#include "pism/icebin/NullTransportHydrology.hh"
#include "pism/geometry/Geometry.hh"
namespace pism {
namespace icebin {


NullTransportHydrology::NullTransportHydrology(std::shared_ptr<const pism::Grid> grid)
    : pism::hydrology::NullTransport(grid), basal_runoff_sum(m_grid, "basal_runoff") {
  basal_runoff_sum.metadata(0)
      .long_name(
          "Effective thickness of subglacial water expelled from till (thickness of water, not ice)")
      .units("m");
}


//! Update the till water thickness by simply integrating the melt input.
/*!
Does an explicit (Euler) step of the integration
  \f[ \frac{\partial W_{til}}{\partial t} = \frac{m}{\rho_w} - C\f]
where \f$C=\f$`hydrology_tillwat_decay_rate_null`.  Enforces bounds
\f$0 \le W_{til} \le W_{til}^{max}\f$ where the upper bound is
`hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.

There is no attempt to report on conservation errors because the model does
not conserve water.

There is no tranportable water thickness variable and no interaction with it.
 */
void NullTransportHydrology::update_impl(double t, double dt,
                                         const hydrology::Inputs &inputs) {
  (void) t;
  const double tillwat_max = m_config->get_number("hydrology.tillwat_max"),
               C           = m_config->get_number("hydrology.tillwat_decay_rate");

  if (tillwat_max < 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "hydrology::NullTransport: hydrology.tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  const auto &cell_type = inputs.geometry->cell_type;

  array::AccessScope scope{ &cell_type, &m_Wtill, &m_surface_input_rate, &basal_runoff_sum };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j) || cell_type.ice_free(i, j)) {
      m_Wtill(i, j) = 0.0;
    } else {
      m_Wtill(i, j) += dt * (m_surface_input_rate(i, j) - C);
      auto Wtil0 = m_Wtill(i, j); // ICEBIN ADDITION
      m_Wtill(i, j) = std::min(std::max(0.0, m_Wtill(i, j)), tillwat_max);
      basal_runoff_sum(i, j) += Wtil0 - m_Wtill(i, j); // ICEBIN ADDITION
    }
  }
}
} // end of namespace icebin
} // end of namespace pism
