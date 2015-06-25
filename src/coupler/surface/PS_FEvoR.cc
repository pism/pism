/* Copyright (C) 2014 PISM Authors
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

#include "PS_FEvoR.hh"
#include "base/util/PISMConfigInterface.hh"

namespace pism {

  namespace surface{
    PS_FEvoR::PS_FEvoR(IceGrid::ConstPtr g)
  : PSFormulas(g) {
  // empty
}

PS_FEvoR::~PS_FEvoR() {
  // empty
}


void PS_FEvoR::init(Vars &vars) {

  (void) vars;

  m_log->message(2, "* Initializing PISM-FEvoR climate inputs...\n"); 

  m_climatic_mass_balance.set(0.0); 

  m_ice_surface_temp.set(243.15);

  // convert from [m/s] to [kg m-2 s-1]
  m_climatic_mass_balance.scale(m_config->get_double("ice_density")); 

}

void PS_FEvoR::update(PetscReal t, PetscReal dt) {
  (void) t;
  (void) dt;

  // do nothing (but an implementation is required)

}
  } // end of namespace surface
} // end of namespace pism
