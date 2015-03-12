// Copyright (C) 2011, 2012, 2013, 2014, 2015 Constantine Khroulev
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
#include "IceGrid.hh"

namespace pism {

ConstantYieldStress::ConstantYieldStress(const IceGrid &g)
  : YieldStress(g) {
  // empty
}

ConstantYieldStress::~ConstantYieldStress () {
  // empty
}

void ConstantYieldStress::init_impl() {
  verbPrintf(2, m_grid.com, "* Initializing the constant basal yield stress model...\n");

  options::String i("-i", "PISM input file"),
    bootstrap("-boot_file", "PISM bootstrapping file");

  options::Real
    tauc("-tauc", "set basal yield stress to a constant (units of Pa)",
         m_config.get_double("default_tauc"));

  // if -tauc was set we just use that value
  if (tauc.is_set()) {
    m_tauc.set(tauc);
  } else if (i.is_set() || bootstrap.is_set()) {
    std::string filename;
    int start;
    bool boot = false;
    find_pism_input(filename, boot, start);

    if (i.is_set()) {
      m_tauc.read(filename, start);
    } else {
      m_tauc.regrid(filename, OPTIONAL,
                  m_config.get_double("default_tauc"));
    }
  } else {
    m_tauc.set(m_config.get_double("default_tauc"));
  }

  regrid("ConstantYieldStress", m_tauc);
}

MaxTimestep ConstantYieldStress::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}


void ConstantYieldStress::add_vars_to_output_impl(const std::string &/*keyword*/,
                                                  std::set<std::string> &result) {
  result.insert("tauc");
}


void ConstantYieldStress::define_variables_impl(const std::set<std::string> &vars,
                                                const PIO &nc, IO_Type nctype) {
  if (set_contains(vars, "tauc")) {
    m_tauc.define(nc, nctype);
  }
}


void ConstantYieldStress::write_variables_impl(const std::set<std::string> &vars,
                                               const PIO &nc) {
  if (set_contains(vars, "tauc")) {
    m_tauc.write(nc);
  }
}


void ConstantYieldStress::update_impl(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
}

} // end of namespace pism
