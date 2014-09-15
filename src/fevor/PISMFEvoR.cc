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

#include "PISMFEvoR.hh"
#include "PISMConfig.hh"
#include "PISMVars.hh"

namespace pism {

PISMFEvoR::PISMFEvoR(IceGrid &g, const Config &conf)
  : Component_TS(g, conf) {
  PetscErrorCode ierr = allocate(); CHKERRCONTINUE(ierr);
}

PISMFEvoR::~PISMFEvoR() {
  // empty
}

PetscErrorCode PISMFEvoR::max_timestep(double t, double &dt, bool &restrict) {
  // FIXME: put real code here
  PetscErrorCode ierr = Component_TS::max_timestep(t, dt, restrict); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMFEvoR::update(double t, double dt) {
  m_t = t;
  m_dt = dt;

  // FIXME: put real code here
  PetscErrorCode ierr = m_enhancement_factor.set(config.get("sia_enhancement_factor"));

  return 0;
}

void PISMFEvoR::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword != "none") {
    result.insert(m_enhancement_factor.metadata().get_string("short_name"));
  }
}

PetscErrorCode PISMFEvoR::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                           IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "enhancement_factor")) {
    ierr = m_enhancement_factor.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMFEvoR::write_variables(const std::set<std::string> &vars, const PIO& nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "enhancement_factor")) {
    ierr = m_enhancement_factor.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMFEvoR::allocate() {
  PetscErrorCode ierr;

  // SIAFD diffusive flux computation requires stencil width of 1
  const unsigned int stencil_width = 1;

  ierr = m_enhancement_factor.create(grid, "enhancement_factor", WITH_GHOSTS,
                                     stencil_width); CHKERRQ(ierr);
  ierr = m_enhancement_factor.set_attrs("diagnostic", // i.e. not needed to re-start the model
                                        "flow law enhancement factor",
                                        "1", // dimensionless
                                        "" /* no standard name */); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMFEvoR::init(Vars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "* Initializing the Fabric Evolution with Recrystallization"
                    " (FEvoR) model...\n"); CHKERRQ(ierr);

  // make enhancement factor available to other PISM components:
  ierr = vars.add(m_enhancement_factor); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
