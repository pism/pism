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
#include "PISMAtmosphere.hh"
#include "PISMConfig.hh"
#include "pism_options.hh"

namespace pism {

PS_FEvoR::PS_FEvoR(IceGrid &g, const Config &conf)
  : PS_EISMINTII(g, conf, (int)'A') {
  // empty
}

PS_FEvoR::~PS_FEvoR() {
  // empty
}


PetscErrorCode PS_FEvoR::init(Vars &vars) {
  PetscErrorCode ierr;

  (void) vars;

  ierr = verbPrintf(2, grid.com, "* Initializing PISM-FEvoR climate inputs...\n"); CHKERRQ(ierr);

  ierr = m_climatic_mass_balance.set(0.0); CHKERRQ(ierr);

  ierr = m_ice_surface_temp.set(243.15); CHKERRQ(ierr);

  // convert from [m/s] to [kg m-2 s-1]
  ierr = m_climatic_mass_balance.scale(config.get("ice_density")); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
