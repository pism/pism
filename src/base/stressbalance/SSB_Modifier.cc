// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev and Ed Bueler
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
#include "flowlaws.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "flowlaw_factory.hh"
#include "PISMConfig.hh"

namespace pism {

PetscErrorCode SSB_Modifier::allocate() {
  PetscErrorCode ierr;

  u.create(grid, "uvel", WITH_GHOSTS);
  u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
              "m s-1", "land_ice_x_velocity");
  u.set_glaciological_units("m year-1");
  u.write_in_glaciological_units = true;

  v.create(grid, "vvel", WITH_GHOSTS);
  v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
              "m s-1", "land_ice_y_velocity");
  v.set_glaciological_units("m year-1");
  v.write_in_glaciological_units = true;

  strain_heating.create(grid, "strainheat", WITHOUT_GHOSTS); // never diff'ed in hor dirs
  strain_heating.set_attrs("internal",
                           "rate of strain heating in ice (dissipation heating)",
                           "W m-3", "");
  strain_heating.set_glaciological_units("mW m-3");

  diffusive_flux.create(grid, "diffusive_flux", WITH_GHOSTS, 1);
  diffusive_flux.set_attrs("internal", 
                           "diffusive (SIA) flux components on the staggered grid",
                           "", "");

  return 0;
}

void ConstantInColumn::init(Vars &vars) {

  SSB_Modifier::init(vars);
}

ConstantInColumn::ConstantInColumn(IceGrid &g, EnthalpyConverter &e, const Config &c)
  : SSB_Modifier(g, e, c)
{
  IceFlowLawFactory ice_factory(grid.com, "sia_", config, &EC);

  ice_factory.setType(config.get_string("sia_flow_law"));

  ice_factory.setFromOptions();
  ice_factory.create(&flow_law);
}

ConstantInColumn::~ConstantInColumn()
{
  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
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
void ConstantInColumn::update(IceModelVec2V *vel_input, bool fast) {

  if (fast) {
    return;
  }

  // horizontal velocity and its maximum:
  IceModelVec::AccessList list;
  list.add(u);
  list.add(v);
  list.add(*vel_input);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    u.setColumn(i,j, (*vel_input)(i,j).u);
    v.setColumn(i,j, (*vel_input)(i,j).v);
  }

  // Communicate to get ghosts (needed to compute w):
  u.update_ghosts();
  v.update_ghosts();

  // diffusive flux and maximum diffusivity
  diffusive_flux.set(0.0);
  D_max = 0.0;
}

} // end of namespace pism
