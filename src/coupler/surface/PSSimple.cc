// Copyright (C) 2011 PISM Authors
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

#include "PSSimple.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "iceModelVec.hh"

///// Simple PISM surface model.

PetscErrorCode PSSimple::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(grid.com, 1, "PISMSurfaceModel::init(PISMVars &vars): atmosphere == NULL");

  ierr = atmosphere->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "* Initializing the simplest PISM surface (snow) processes model PSSimple.\n"
     "  It passes atmospheric state directly to upper ice fluid surface:\n"
     "    surface mass balance          := precipitation,\n"
     "    ice upper surface temperature := 2m air temperature.\n");
     CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSSimple::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_precip(result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted precipitation as surface mass balance (PSSimple)\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean annual 2 m air temperature as instantaneous ice temperature at the ice surface (PSSimple)\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

void PSSimple::add_vars_to_output(string keyword, set<string> &result) {
  atmosphere->add_vars_to_output(keyword, result);
}
