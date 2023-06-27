// Copyright (C) 2012-2014, 2016 PISM Authors
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

#include <base/hydrology/hydrology_diagnostics.hh>
#include <base/util/IceModelVec2CellType.hh>
#include <base/util/Mask.hh>
#include <base/util/PISMVars.hh>
#include <base/util/error_handling.hh>
#include <icebin/NullTransportHydrology.hh>

namespace pism {
namespace icebin {


NullTransportHydrology::NullTransportHydrology(pism::IceGrid::ConstPtr g) : pism::hydrology::NullTransport(g) {

  printf("BEGIN NullTransportHydrology::NullTransportHydrology()\n");

  // *all* PISMHydrology classes have layer of water stored in till
  basal_runoff_sum.create(m_grid, "basal_runoff", WITHOUT_GHOSTS);
  basal_runoff_sum.set_attrs("excess water", "Effective thickness of subglacial water expelled from till "
                                             "(thickness of water, not ice)",
                             "m s-1", "");

  printf("END NullTransportHydrology::NullTransportHydrology()\n");
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
void NullTransportHydrology::update_impl(double icet, double icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - m_t) < 1e-6) && (fabs(icedt - m_dt) < 1e-6)) {
    return;
  }
  m_t  = icet;
  m_dt = icedt;

  get_input_rate(icet, icedt, m_total_input);

  const double tillwat_max = m_config->get_double("hydrology.tillwat_max"),
               C           = m_config->get_double("hydrology.tillwat_decay_rate");

  if (tillwat_max < 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "hydrology::NullTransport: hydrology.tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  const IceModelVec2CellType &cell_type = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(m_Wtil);
  list.add(m_total_input);
  list.add(basal_runoff_sum); // ICEBIN ADDITION
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j) || cell_type.ice_free(i, j)) {
      m_Wtil(i, j) = 0.0;
    } else {
      m_Wtil(i, j) += icedt * (m_total_input(i, j) - C);
      auto Wtil0 = m_Wtil(i, j); // ICEBIN ADDITION
      m_Wtil(i, j) = std::min(std::max(0.0, m_Wtil(i, j)), tillwat_max);
      basal_runoff_sum(i, j) += Wtil0 - m_Wtil(i, j); // ICEBIN ADDITION
    }
  }
}
} // end of namespace icebin
} // end of namespace pism
