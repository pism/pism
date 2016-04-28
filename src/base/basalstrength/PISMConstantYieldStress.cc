// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#include "base/util/pism_options.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/IceGrid.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

ConstantYieldStress::ConstantYieldStress(IceGrid::ConstPtr g)
  : YieldStress(g) {
  // empty
}

ConstantYieldStress::~ConstantYieldStress () {
  // empty
}

void ConstantYieldStress::init_impl() {
  m_log->message(2, "* Initializing the constant basal yield stress model...\n");

  std::string filename;
  int start = 0;
  bool boot = false;
  bool use_input_file = find_pism_input(filename, boot, start);
  double tauc = m_config->get_double("basal_yield_stress.constant.value");
  if (use_input_file) {
    if (boot) {
      m_tauc.regrid(filename, OPTIONAL, tauc);
    } else {
      m_tauc.read(filename, start);
    }
  } else {
    // Set the constant value.
    m_tauc.set(tauc);
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
