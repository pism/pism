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

#include "PSLapseRates.hh"

PetscErrorCode PSLapseRates::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool smb_lapse_rate_set;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "  [using temperature and mass balance lapse corrections]\n"); CHKERRQ(ierr);

  ierr = init_internal(vars); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Lapse rate options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsReal("-smb_lapse_rate",
                           "Elevation lapse rate for the surface mass balance, in m/year per km",
                           smb_lapse_rate, smb_lapse_rate_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "   ice upper-surface temperature lapse rate: %3.3f K per km\n"
                    "   ice-equivalent surface mass balance lapse rate: %3.3f m/year per km\n",
                    temp_lapse_rate, smb_lapse_rate); CHKERRQ(ierr);

  temp_lapse_rate = convert(temp_lapse_rate, "K/km", "K/m");

  smb_lapse_rate = convert(smb_lapse_rate, "m/year / km", "m/s / m");

  return 0;
}

PetscErrorCode PSLapseRates::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->ice_surface_mass_flux(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, smb_lapse_rate); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSLapseRates::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, temp_lapse_rate); CHKERRQ(ierr);
  return 0;
}
