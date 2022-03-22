// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2022 Constantine Khroulev and Ed Bueler
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

#include "SSB_Modifier.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Vars.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace stressbalance {

SSB_Modifier::SSB_Modifier(IceGrid::ConstPtr g)
  : Component(g),
    m_EC(g->ctx()->enthalpy_converter()),
    m_diffusive_flux(m_grid, "diffusive_flux"),
    m_u(m_grid, "uvel", array::WITH_GHOSTS, m_grid->z()),
    m_v(m_grid, "vvel", array::WITH_GHOSTS, m_grid->z()) {
  m_D_max = 0.0;

  m_u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
                "m s-1", "m year-1", "land_ice_x_velocity", 0);

  m_v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
                "m s-1", "m year-1", "land_ice_y_velocity", 0);

  m_diffusive_flux.set_attrs("internal",
                             "diffusive (SIA) flux components on the staggered grid",
                             "", "", "", 0);

}

void SSB_Modifier::init() {
}

const array::Staggered& SSB_Modifier::diffusive_flux() {
  return m_diffusive_flux;
}

//! \brief Get the max diffusivity (for the adaptive time-stepping).
double SSB_Modifier::max_diffusivity() const {
  return m_D_max;
}

const array::Array3D& SSB_Modifier::velocity_u() const {
  return m_u;
}

const array::Array3D& SSB_Modifier::velocity_v() const {
  return m_v;
}

std::string SSB_Modifier::stdout_report() const {
  return "";
}

std::shared_ptr<const rheology::FlowLaw> SSB_Modifier::flow_law() const {
  return m_flow_law;
}

void ConstantInColumn::init() {
  SSB_Modifier::init();
}

ConstantInColumn::ConstantInColumn(IceGrid::ConstPtr g)
  : SSB_Modifier(g) {
  rheology::FlowLawFactory ice_factory("stress_balance.sia.", m_config, m_EC);

  m_flow_law = ice_factory.create();
}


//! \brief Distribute the input velocity throughout the column.
/*!
 * Things to update:
 * - 3D-distributed horizontal velocity
 * - maximum horizontal velocity
 * - diffusive ice flux
 * - maximum diffusivity
 * - strain heating (strain_heating)
 */
void ConstantInColumn::update(const array::Vector &sliding_velocity,
                              const Inputs &inputs,
                              bool full_update) {

  (void) inputs;

  if (not full_update) {
    return;
  }

  // horizontal velocity and its maximum:
  array::AccessScope list{&m_u, &m_v, &sliding_velocity};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_u.set_column(i,j, sliding_velocity(i,j).u);
    m_v.set_column(i,j, sliding_velocity(i,j).v);
  }

  // Communicate to get ghosts (needed to compute w):
  m_u.update_ghosts();
  m_v.update_ghosts();

  // diffusive flux and maximum diffusivity
  m_diffusive_flux.set(0.0);
  m_D_max = 0.0;
}

} // end of namespace stressbalance
} // end of namespace pism
