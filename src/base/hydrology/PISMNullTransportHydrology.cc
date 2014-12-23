// Copyright (C) 2012-2014 PISM Authors
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

#include "Mask.hh"
#include "PISMHydrology.hh"
#include "hydrology_diagnostics.hh"
#include "error_handling.hh"

namespace pism {

NullTransportHydrology::NullTransportHydrology(const IceGrid &g)
  : Hydrology(g) {
}

NullTransportHydrology::~NullTransportHydrology() {
}

void NullTransportHydrology::init() {
  verbPrintf(2, m_grid.com,
             "* Initializing the null-transport (till only) subglacial hydrology model ...\n");
  Hydrology::init();
}


//! Set the transportable subglacial water thickness to zero; there is no tranport.
void NullTransportHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  result.set(0.0);
}


//! Returns the (trivial) overburden pressure as the pressure of the non-existent transportable water, because this is the least harmful output if this is misused.
void NullTransportHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  overburden_pressure(result);
}


//! Update the till water thickness by simply integrating the melt input.
/*!
Does a step of the trivial integration
  \f[ \frac{\partial W_{til}}{\partial t} = \frac{m}{\rho_w} - C\f]
where \f$C=\f$`hydrology_tillwat_decay_rate`.  Enforces bounds
\f$0 \le W_{til} \le W_{til}^{max}\f$ where the upper bound is
`hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.

Uses the current mass-continuity timestep `icedt`.  (Compare
RoutingHydrology::raw_update_Wtil() which will generally be taking time steps
determined by the evolving transportable water layer in that model.)

There is no attempt to report on conservation errors because this
NullTransportHydrology model does not conserve water.

There is no tranportable water thickness variable and no interaction with it.
 */
void NullTransportHydrology::update(double icet, double icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - m_t) < 1e-6) && (fabs(icedt - m_dt) < 1e-6)) {
    return;
  }
  m_t = icet;
  m_dt = icedt;

  get_input_rate(icet,icedt,total_input);

  const double tillwat_max = m_config.get("hydrology_tillwat_max"),
               C           = m_config.get("hydrology_tillwat_decay_rate");

  if (tillwat_max < 0.0) {
    throw RuntimeError("NullTransportHydrology: hydrology_tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  MaskQuery M(*mask);
  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(Wtil);
  list.add(total_input);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ocean(i,j) || M.ice_free(i,j)) {
      Wtil(i,j) = 0.0;
    } else {
      Wtil(i,j) += icedt * (total_input(i,j) - C);
      Wtil(i,j) = std::min(std::max(0.0, Wtil(i,j)), tillwat_max);
    }
  }
}


} // end of namespace pism
