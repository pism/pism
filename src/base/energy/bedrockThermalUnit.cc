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


PISMBedThermalUnit::PISMBedThermalUnit(IceGrid &g, EnthalpyConverter &e, 
                                       const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf), EC(e) {
  mask = NULL;
  thk = NULL;
  enthalpy = NULL;

  if (allocate()) {
    verbPrintf(1,g.com, "PISMBedThermalUnit::allocate() returned nonzero\n");
    PISMEnd();
  }
}

PetscErrorCode PISMBedThermalUnit::allocate() {
  PetscErrorCode ierr;

  // FIXME:  this temp object will *not* be creatable from IceGrid once
  //         the PISMBedThermalUnit is seriously implemented; instead it will read its own
  //         -Mbz (and -Lbz) and create a dof=Mbz DA accordingly, presumably
  ierr = temp.create(grid, "btu_litho_temp", false); CHKERRQ(ierr);
  ierr = temp.set_attrs("model_state",
                        "lithosphere (bedrock) temperature, in PISMBedThermalUnit",
		        "K", ""); CHKERRQ(ierr);
  ierr = temp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  
  ierr = ghf.create(grid, "btu_bheatflx", false); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  ierr = ghf.set_attrs("climate_steady",
                       "upward geothermal flux at bedrock thermal layer base, deep in lithosphere",
		       "W m-2", ""); CHKERRQ(ierr);
  ierr = ghf.set_glaciological_units("mW m-2");
  ghf.time_independent = true;

  ierr = ice_base_temp.create(grid, "btu_ice_base_temp", false); CHKERRQ(ierr);
  ierr = ice_base_temp.set_attrs("internal",
    "temperature at base of ice for duration of timestep, in PISMBedThermalUnit",
    "K", ""); CHKERRQ(ierr);
  ierr = ice_base_temp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode PISMBedThermalUnit::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com,"* Initializing the bedrock thermal unit ...\n"); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");
  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(2, "thk is not available");
  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(3, "enthalpy is not available");

  ice_base_temp.set(-1.0); // init as invalid

  // build constant diffusivity for heat equation
  bed_rho = config.get("bedrock_thermal_density");
  bed_c   = config.get("bedrock_thermal_specific_heat_capacity");
  bed_k   = config.get("bedrock_thermal_conductivity");
  bed_D   = bed_k / (bed_rho * bed_c);

  // read state information from file, or regrid it
  bool regrid = false;
  int start = -1;
  string filename;
  ierr = find_pism_input(filename, regrid, start); CHKERRQ(ierr);
  string tempname = temp.string_attr("short_name"),
         ghfname = ghf.string_attr("short_name");
  ierr = verbPrintf(2, grid.com, 
                    "    reading %s and %s from file %s ... \n",
                    tempname.c_str(),ghfname.c_str(),filename.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = temp.regrid(filename.c_str(), true); CHKERRQ(ierr); // fails if not found!
    ierr = ghf.regrid(filename.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = temp.read(filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = ghf.read(filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }
  string history = "read from " + filename + "\n";
  ierr = temp.set_attr("history", history); CHKERRQ(ierr);
  ierr = ghf.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}


// FIXME:  bootstrapping code is needed when running with IceModel, at least


void PISMBedThermalUnit::add_vars_to_output(string keyword, set<string> &result) {
  result.insert(temp.string_attr("short_name"));
  result.insert(ghf.string_attr("short_name"));
  if (keyword == "big") {
    result.insert(ice_base_temp.string_attr("short_name"));
  }
}


PetscErrorCode PISMBedThermalUnit::define_variables(
                         set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, temp.string_attr("short_name"))) {
    ierr = temp.define(nc, nctype); CHKERRQ(ierr); 
  }  
  if (set_contains(vars, ghf.string_attr("short_name"))) {
    ierr = ghf.define(nc, nctype); CHKERRQ(ierr); 
  }  
  if (set_contains(vars, ice_base_temp.string_attr("short_name"))) {
    ierr = ice_base_temp.define(nc, nctype); CHKERRQ(ierr); 
  }  
  return 0;
}


PetscErrorCode PISMBedThermalUnit::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  if (set_contains(vars, temp.string_attr("short_name"))) {
    ierr = temp.write(filename.c_str()); CHKERRQ(ierr); 
  }  
  if (set_contains(vars, ghf.string_attr("short_name"))) {
    ierr = ghf.write(filename.c_str()); CHKERRQ(ierr); 
  }  
  if (set_contains(vars, ice_base_temp.string_attr("short_name"))) {
    ierr = ice_base_temp.write(filename.c_str()); CHKERRQ(ierr); 
  }  
  return 0;
}


PetscErrorCode PISMBedThermalUnit::max_timestep(PetscReal /*t_years*/, PetscReal &dt_years) {

  const PetscReal dzb = grid.dz_fine;
  dt_years = dzb * dzb / (2.0 * bed_D);
  dt_years /= secpera;
  return 0;
}


PetscErrorCode PISMBedThermalUnit::get_upward_geothermal_flux(
                          PetscReal t_years, PetscReal dt_years,
                          IceModelVec2S &result) {
  PetscErrorCode ierr;

  PetscScalar mydtyears;
  ierr = max_timestep(t_years,mydtyears); CHKERRQ(ierr);
  if (mydtyears < dt_years) {
     SETERRQ(4,"PISMBedThermalUnit::get_upward_geothermal_flux() thinks you asked\n"
               "    for too big a timestep\n");
  }

  //if (mask == NULL) { SETERRQ(1, "mask unavailable"); }
  if (thk == NULL) { SETERRQ(2, "thk unavailable"); }
  if (enthalpy == NULL) { SETERRQ(3, "enthalpy unavailable"); }

  // first job: fill ice_base_temp; FIXME: is correct in floating case? ice-free case?
  ierr = enthalpy->getHorSlice(ice_base_temp, 0.0); CHKERRQ(ierr);
  ierr = ice_base_temp.begin_access(); CHKERRQ(ierr);
  ierr = thk->begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      const PetscReal pressure = EC.getPressureFromDepth((*thk)(i,j));
      ierr = EC.getAbsTemp(ice_base_temp(i,j), pressure, ice_base_temp(i,j));
               CHKERRQ(ierr);
    }
  }
  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = ice_base_temp.end_access(); CHKERRQ(ierr);
  // no need to communicate just-filled ice_base_temp

  const PetscReal dzb = grid.dz_fine;
  const PetscInt  Mbz = grid.Mbz_fine,
                  k0  = Mbz - 1;          // Tb[k0] = ice/bed interface temp, at z=0

  const PetscReal bed_R  = bed_D * (dt_years * secpera) / (dzb * dzb);

  if (Mbz <= 1) {
    SETERRQ(5,"Mbz <= 1 not allowed in PISMBedThermalUnit::get_upward_geothermal_flux()\n"); }

  PetscScalar *Tbold, *Tbnew;
  Tbold = new PetscScalar[Mbz];
  Tbnew = new PetscScalar[Mbz];

  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = ghf.begin_access(); CHKERRQ(ierr);
  ierr = ice_base_temp.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      ierr = temp.getValColumnPL(i,j,Tbold); CHKERRQ(ierr);
      Tbold[k0] = ice_base_temp(i,j);  // sets Dirichlet explicit-in-time at top of column

      const PetscReal Tbold_negone = Tbold[1] + 2 * ghf(i,j) * dzb / bed_k;
      Tbnew[0] = Tbold[0] + bed_R * (Tbold_negone - 2 * Tbold[0] + Tbold[1]);
      for (PetscInt k = 1; k < k0; k++) { // working upward from base
        Tbnew[k] = Tbold[k] + bed_R * (Tbold[k-1] - 2 * Tbold[k] + Tbold[k+1]);
      }
      Tbnew[k0] = ice_base_temp(i,j);

      ierr = temp.setValColumnPL(i,j,Tbnew); CHKERRQ(ierr);

      result(i,j) = - bed_k * (Tbnew[k0] - Tbnew[k0-1]) / dzb;  
    }
  }
  ierr = ice_base_temp.end_access(); CHKERRQ(ierr);
  ierr = ghf.end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  delete [] Tbold;
  delete [] Tbnew;

  // FIXME: a derived class version of this method will be simpler: will just report the values of the stored geothermal flux

  return 0;
}

