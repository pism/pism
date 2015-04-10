/* Copyright (C) 2014, 2015 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "Context.hh"

namespace pism {

class Context::Impl {
public:
  Impl(MPI_Comm c,
       UnitSystemPtr sys,
       ConfigPtr conf,
       EnthalpyConverterPtr EC,
       TimePtr t)
    : com(c), unit_system(sys), config(conf), enthalpy_converter(EC), time(t) {
    // empty
  }
  MPI_Comm com;
  UnitSystemPtr unit_system;
  ConfigPtr config;
  EnthalpyConverterPtr enthalpy_converter;
  TimePtr time;
};

Context::Context(MPI_Comm c, UnitSystemPtr sys,
                 ConfigPtr conf, EnthalpyConverterPtr EC, TimePtr t)
  : m_impl(new Impl(c, sys, conf, EC, t)) {
  // empty
}

MPI_Comm Context::com() const {
  return m_impl->com;      
}

Context::UnitSystemPtr Context::unit_system() const {
  return m_impl->unit_system;
}

Context::ConfigPtr Context::config() {
  return m_impl->config;
}

Context::ConstConfigPtr Context::config() const {
  return m_impl->config;
}

Context::EnthalpyConverterPtr Context::enthalpy_converter() const {
  return m_impl->enthalpy_converter;
}

Context::TimePtr Context::time() {
  return m_impl->time;     
}

Context::ConstTimePtr Context::time() const {
  return m_impl->time;    
}

} // end of namespace pism
