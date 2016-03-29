// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev and Ed Bueler
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
#include "base/rheology/FlowLawFactory.hh"
#include "base/rheology/FlowLaw.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMVars.hh"

namespace pism {
namespace stressbalance {

SSB_Modifier::SSB_Modifier(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e)
  : Component(g), m_EC(e) {

  m_D_max = 0.0;

  m_u.create(m_grid, "uvel", WITH_GHOSTS);
  m_u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
              "m s-1", "land_ice_x_velocity");
  m_u.metadata().set_string("glaciological_units", "m year-1");
  m_u.write_in_glaciological_units = true;

  m_v.create(m_grid, "vvel", WITH_GHOSTS);
  m_v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
              "m s-1", "land_ice_y_velocity");
  m_v.metadata().set_string("glaciological_units", "m year-1");
  m_v.write_in_glaciological_units = true;

  m_strain_heating.create(m_grid, "strainheat", WITHOUT_GHOSTS); // never diff'ed in hor dirs
  m_strain_heating.set_attrs("internal",
                           "rate of strain heating in ice (dissipation heating)",
                           "W m-3", "");
  m_strain_heating.metadata().set_string("glaciological_units", "mW m-3");

  m_diffusive_flux.create(m_grid, "diffusive_flux", WITH_GHOSTS, 1);
  m_diffusive_flux.set_attrs("internal", 
                           "diffusive (SIA) flux components on the staggered grid",
                           "", "");
  
}

SSB_Modifier::~SSB_Modifier() {
  // empty
}

void SSB_Modifier::init() {
}

const IceModelVec2Stag& SSB_Modifier::diffusive_flux() {
  return m_diffusive_flux;
}

//! \brief Get the max diffusivity (for the adaptive time-stepping).
double SSB_Modifier::max_diffusivity() {
  return m_D_max;
}

const IceModelVec3& SSB_Modifier::velocity_u() {
  return m_u;
}

const IceModelVec3& SSB_Modifier::velocity_v() {
  return m_v;
}

const IceModelVec3& SSB_Modifier::volumetric_strain_heating() {
  return m_strain_heating;
}

std::string SSB_Modifier::stdout_report() {
  return "";
}

rheology::FlowLaw* SSB_Modifier::flow_law() {
  return m_flow_law;
}


void ConstantInColumn::init() {
  SSB_Modifier::init();
}

ConstantInColumn::ConstantInColumn(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e)
  : SSB_Modifier(g, e)
{
  rheology::FlowLawFactory ice_factory("sia_", m_config, m_EC);

  m_flow_law = ice_factory.create();
}

ConstantInColumn::~ConstantInColumn()
{
  if (m_flow_law != NULL) {
    delete m_flow_law;
    m_flow_law = NULL;
  }
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
void ConstantInColumn::update(const IceModelVec2V &vel_input, bool fast) {

  if (fast) {
    return;
  }

  // horizontal velocity and its maximum:
  IceModelVec::AccessList list;
  list.add(m_u);
  list.add(m_v);
  list.add(vel_input);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_u.set_column(i,j, vel_input(i,j).u);
    m_v.set_column(i,j, vel_input(i,j).v);
  }

  // Communicate to get ghosts (needed to compute w):
  m_u.update_ghosts();
  m_v.update_ghosts();

  // diffusive flux and maximum diffusivity
  m_diffusive_flux.set(0.0);
  m_D_max = 0.0;
}

void ConstantInColumn::add_vars_to_output_impl(const std::string &keyword,
                                          std::set<std::string> &result) {
  (void)keyword;
  (void)result;
}

void ConstantInColumn::define_variables_impl(const std::set<std::string> &vars,
                                        const PIO &nc,
                                        IO_Type nctype) {
  (void)vars;
  (void)nc;
  (void)nctype;
}

void ConstantInColumn::write_variables_impl(const std::set<std::string> &vars,
                                       const PIO &nc) {
  (void)vars;
  (void)nc;
}


} // end of namespace stressbalance
} // end of namespace pism
