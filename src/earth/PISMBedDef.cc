// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
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

#include "PISMBedDef.hh"
#include "PISMTime.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"

namespace pism {

BedDef::BedDef(IceGrid &g, const Config &conf)
  : Component_TS(g, conf) {

  thk    = NULL;
  topg   = NULL;
  uplift = NULL;

  t_beddef_last = GSL_NAN;

  PetscErrorCode ierr = pismbeddef_allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "BedDef::BedDef(...): pismbeddef_allocate() failed\n");
    PISMEnd();
  }
}

PetscErrorCode BedDef::pismbeddef_allocate() {
  PetscErrorCode ierr;
  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");

  ierr = topg_initial.create(grid, "topg_initial", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_initial.set_attrs("model_state", "bedrock surface elevation (at the beginning of the run)",
                                "m", ""); CHKERRQ(ierr);

  ierr = topg_last.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_last.set_attrs("model_state", "bedrock surface elevation",
                             "m", "bedrock_altitude"); CHKERRQ(ierr);

  return 0;
}

void BedDef::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("topg_initial");
}

PetscErrorCode BedDef::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                            IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "topg_initial")) {
    ierr = topg_initial.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode BedDef::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "topg_initial")) {
    ierr = topg_initial.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode BedDef::init(Vars &vars) {
  PetscErrorCode ierr;

  t_beddef_last = grid.time->start();

  thk    = vars.get_2d_scalar("land_ice_thickness");
  topg   = vars.get_2d_scalar("bedrock_altitude");
  uplift = vars.get_2d_scalar("tendency_of_bedrock_altitude");

  // Save the bed elevation at the beginning of the run:
  ierr = topg_initial.copy_from(*topg); CHKERRQ(ierr);

  return 0;
}

//! Compute bed uplift (dt_beddef is in seconds).
PetscErrorCode BedDef::compute_uplift(double dt_beddef) {
  PetscErrorCode ierr;

  ierr = topg->add(-1, topg_last, *uplift); CHKERRQ(ierr);
  //! uplift = (topg - topg_last) / dt
  ierr = uplift->scale(1.0 / dt_beddef); CHKERRQ(ierr); 

  return 0;
}

} // end of namespace pism
