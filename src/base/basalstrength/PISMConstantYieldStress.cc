// Copyright (C) 2011, 2012 Constantine Khroulev
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

#include "PISMConstantYieldStress.hh"
#include "pism_options.hh"

PetscErrorCode PISMConstantYieldStress::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  bool i_set, bootstrap, tauc_set;
  PetscReal constant_tauc;
  string filename;
  int start;

  ierr = verbPrintf(2, grid.com, "* Initializing the constant basal yield stress model...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "PISMConstantYieldStress options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-tauc", "set basal yield stress to a constant (units of Pa)",
                           constant_tauc, tauc_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // if -tauc was set we just use that value
  if (tauc_set) {
    ierr = tauc.set(constant_tauc); CHKERRQ(ierr);
  } else if (i_set || bootstrap) {
    ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);

    if (i_set) {
      ierr = tauc.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = tauc.regrid(filename, config.get("default_tauc")); CHKERRQ(ierr);
    }
  } else {
    ierr = tauc.set(config.get("default_tauc")); CHKERRQ(ierr);
  }

  ierr = regrid(); CHKERRQ(ierr);

  return 0;
}


void PISMConstantYieldStress::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["tauc"] = tauc.get_metadata();
}


PetscErrorCode PISMConstantYieldStress::define_variables(set<string> vars, const PIO &nc,
                                                         PISM_IO_Type nctype) {
  if (set_contains(vars, "tauc")) {
    PetscErrorCode ierr = tauc.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMConstantYieldStress::write_variables(set<string> vars, string filename) {
  if (set_contains(vars, "tauc")) {
    PetscErrorCode ierr = tauc.write(filename); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMConstantYieldStress::update(PetscReal my_t, PetscReal my_dt) {
  t = my_t; dt = my_dt;
  return 0;
}


PetscErrorCode PISMConstantYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  PetscErrorCode ierr = tauc.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMConstantYieldStress::allocate() {
  PetscErrorCode ierr;

  ierr = tauc.create(grid, "tauc", true, grid.max_stencil_width); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  ierr = tauc.set_attrs("model_state", 
                        "yield stress for basal till (plastic or pseudo-plastic model)",
                        "Pa", ""); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMConstantYieldStress::regrid() {
  PetscErrorCode ierr;
  bool regrid_file_set, regrid_vars_set;
  string regrid_file;
  vector<string> regrid_vars;

  ierr = PetscOptionsBegin(grid.com, "", "PISMMohrCoulombYieldStress regridding options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "regridding file name",
                             regrid_file, regrid_file_set); CHKERRQ(ierr);
    ierr = PISMOptionsStringArray("-regrid_vars", "comma-separated list of regridding variables",
                                  "", regrid_vars, regrid_vars_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (! regrid_file_set) return 0;

  set<string> vars;
  for (unsigned int i = 0; i < regrid_vars.size(); ++i)
    vars.insert(regrid_vars[i]);

  // stop if the user did not ask to regrid tillphi
  if (!set_contains(vars, tauc.string_attr("short_name")))
    return 0;

  ierr = tauc.regrid(regrid_file, true); CHKERRQ(ierr);

  return 0;
}

