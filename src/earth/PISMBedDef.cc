// Copyright (C) 2010, 2011, 2012 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

PISMBedDef::PISMBedDef(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent_TS(g, conf) {

  thk    = NULL;
  topg   = NULL;
  uplift = NULL;

  t_beddef_last = GSL_NAN;

  PetscErrorCode ierr = pismbeddef_allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISMBedDef::PISMBedDef(...): pismbeddef_allocate() failed\n");
    PISMEnd();
  }
}

PetscErrorCode PISMBedDef::pismbeddef_allocate() {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = grid.max_stencil_width;

  ierr = topg_initial.create(grid, "topg_initial", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_initial.set_attrs("model_state", "bedrock surface elevation (at the beginning of the run)",
                                "m", ""); CHKERRQ(ierr);

  ierr = topg_last.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_last.set_attrs("model_state", "bedrock surface elevation",
			     "m", "bedrock_altitude"); CHKERRQ(ierr);

  return 0;
}

void PISMBedDef::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["topg_initial"] = topg_initial.get_metadata();
}

PetscErrorCode PISMBedDef::define_variables(set<string> vars, const PIO &nc,
                                            PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "topg_initial")) {
    ierr = topg_initial.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMBedDef::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "topg_initial")) {
    ierr = topg_initial.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMBedDef::init(PISMVars &vars) {
  PetscErrorCode ierr;

  t_beddef_last = grid.time->start();

  thk = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!thk) SETERRQ(grid.com, 1, "ERROR: thk is not available");

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (!topg) SETERRQ(grid.com, 1, "ERROR: topg is not available");

  uplift = dynamic_cast<IceModelVec2S*>(vars.get("tendency_of_bedrock_altitude"));
  if (!uplift) SETERRQ(grid.com, 1, "ERROR: uplift is not available");

  // Save the bed elevation at the beginning of the run:
  ierr = topg_initial.copy_from(*topg); CHKERRQ(ierr);

  return 0;
}

//! Compute bed uplift (dt_beddef is in seconds).
PetscErrorCode PISMBedDef::compute_uplift(PetscScalar dt_beddef) {
  PetscErrorCode ierr;

  ierr = topg->add(-1, topg_last, *uplift); CHKERRQ(ierr);
  //! uplift = (topg - topg_last) / dt
  ierr = uplift->scale(1.0 / dt_beddef); CHKERRQ(ierr); 

  return 0;
}
