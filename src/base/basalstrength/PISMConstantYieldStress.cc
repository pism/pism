// Copyright (C) 2011, 2012, 2013, 2014 Constantine Khroulev
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

#include "PISMConstantYieldStress.hh"
#include "pism_options.hh"
#include "PISMConfig.hh"

PetscErrorCode PISMConstantYieldStress::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  bool i_set, bootstrap, tauc_set;
  double constant_tauc;
  std::string filename;
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
      ierr = tauc.regrid(filename, OPTIONAL,
                         config.get("default_tauc")); CHKERRQ(ierr);
    }
  } else {
    ierr = tauc.set(config.get("default_tauc")); CHKERRQ(ierr);
  }

  ierr = regrid("PISMConstantYieldStress", &tauc); CHKERRQ(ierr);

  return 0;
}


void PISMConstantYieldStress::add_vars_to_output(std::string /*keyword*/, std::set<std::string> &result) {
  result.insert("tauc");
}


PetscErrorCode PISMConstantYieldStress::define_variables(std::set<std::string> vars, const PIO &nc,
                                                         PISM_IO_Type nctype) {
  if (set_contains(vars, "tauc")) {
    PetscErrorCode ierr = tauc.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMConstantYieldStress::write_variables(std::set<std::string> vars, const PIO &nc) {
  if (set_contains(vars, "tauc")) {
    PetscErrorCode ierr = tauc.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMConstantYieldStress::update(double my_t, double my_dt) {
  m_t = my_t; m_dt = my_dt;
  return 0;
}


PetscErrorCode PISMConstantYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  PetscErrorCode ierr = tauc.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMConstantYieldStress::allocate() {
  PetscErrorCode ierr;

  ierr = tauc.create(grid, "tauc", WITH_GHOSTS, grid.max_stencil_width); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  ierr = tauc.set_attrs("model_state", 
                        "yield stress for basal till (plastic or pseudo-plastic model)",
                        "Pa", ""); CHKERRQ(ierr);
  return 0;
}

