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

namespace pism {

void ConstantYieldStress::init(Vars &/*vars*/) {
  bool i_set, bootstrap, tauc_set;
  double constant_tauc = config.get("default_tauc");
  std::string filename;
  int start;

  verbPrintf(2, grid.com, "* Initializing the constant basal yield stress model...\n");

  {
    OptionsIsSet("-i", "PISM input file", i_set);
    OptionsIsSet("-boot_file", "PISM bootstrapping file",
                 bootstrap);
    OptionsReal("-tauc", "set basal yield stress to a constant (units of Pa)",
                constant_tauc, tauc_set);
  }

  // if -tauc was set we just use that value
  if (tauc_set) {
    tauc.set(constant_tauc);
  } else if (i_set || bootstrap) {
    find_pism_input(filename, bootstrap, start);

    if (i_set) {
      tauc.read(filename, start);
    } else {
      tauc.regrid(filename, OPTIONAL,
                  config.get("default_tauc"));
    }
  } else {
    tauc.set(config.get("default_tauc"));
  }

  regrid("ConstantYieldStress", &tauc);
}


void ConstantYieldStress::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("tauc");
}


void ConstantYieldStress::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                                         IO_Type nctype) {
  if (set_contains(vars, "tauc")) {
    tauc.define(nc, nctype);
  }
}


void ConstantYieldStress::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "tauc")) {
    tauc.write(nc);
  }
}


void ConstantYieldStress::update(double my_t, double my_dt) {
  m_t = my_t; m_dt = my_dt;
}


void ConstantYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  tauc.copy_to(result);
}

PetscErrorCode ConstantYieldStress::allocate() {
  tauc.create(grid, "tauc", WITH_GHOSTS, config.get("grid_max_stencil_width"));
  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  tauc.set_attrs("model_state", 
                 "yield stress for basal till (plastic or pseudo-plastic model)",
                 "Pa", "");
  return 0;
}


} // end of namespace pism
