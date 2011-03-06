// Copyright (C) 2011 Ed Bueler
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

#include "bedrockThermalUnit.hh"


BedrockThermalUnit::BedrockThermalUnit(IceGrid &g, EnthalpyConverter &e, 
                                       const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf), EC(e) {
  mask = NULL;
  enthalpy = NULL;

  allocate();
}


PetscErrorCode BedrockThermalUnit::allocate() {
  PetscErrorCode ierr;

  ierr = basal_temp.create(grid, "btu_basal_temp", false); CHKERRQ(ierr);
  ierr = basal_temp.set_attrs("internal",
    "temperature at base of ice for duration of timestep, as presented to the BedrockThermalUnit",
    "K", ""); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode BedrockThermalUnit::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com,"* Initializing the bedrock thermal unit ...\n"); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(1, "enthalpy is not available");

  basal_temp.set(0.0);

  rho = config.get("bedrock_thermal_density");
  c   = config.get("bedrock_thermal_specific_heat_capacity");
  k   = config.get("bedrock_thermal_conductivity");
  
  D = k / (rho * c);

  return 0;
}


PetscErrorCode BedrockThermalUnit::max_timestep(PetscReal /*t_years*/, PetscReal &dt_years) {

  dt_years = grid.dzbMIN * grid.dzbMIN / (2.0 * D);
  dt_years *= secpera;
  return 0;
}


PetscErrorCode BedrockThermalUnit::update(PetscReal /*t_years*/, PetscReal /*dt_years*/) {
  //PetscErrorCode ierr;

  // FIXME: I want to use  IceModel_tempbase  to fill basal_temp here, I believe
  
  // FIXME: then I want to go forward dt_years and update a stored map of 
  //   geothermal flux
  
  // FIXME: note that a derived class will be simpler: this update() will do nothing
  //   and some get_() method will just report the values of the geothermal flux

  SETERRQ(1,"not implemented");
  return 0;
}


PetscErrorCode BedrockThermalUnit::write_model_state(
       PetscReal /*t_years*/, PetscReal /*dt_years*/, string filename) {
  PetscErrorCode ierr;

  ierr = basal_temp.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

