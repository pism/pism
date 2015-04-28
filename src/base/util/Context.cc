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
#include "Profiling.hh"
#include "PISMUnits.hh"
#include "PISMConfig.hh"
#include "PISMTime.hh"
#include "base/enthalpyConverter.hh"

namespace pism {

class Context::Impl {
public:
  Impl(MPI_Comm c,
       UnitsSystemPtr sys,
       ConfigPtr conf,
       EnthalpyConverterPtr EC,
       TimePtr t,
       const std::string &p)
    : com(c), unit_system(sys), config(conf), enthalpy_converter(EC), time(t), prefix(p) {
    // empty
  }
  MPI_Comm com;
  UnitsSystemPtr unit_system;
  ConfigPtr config;
  EnthalpyConverterPtr enthalpy_converter;
  TimePtr time;
  std::string prefix;
  Profiling profiling;
};

Context::Context(MPI_Comm c, UnitsSystemPtr sys,
                 ConfigPtr conf, EnthalpyConverterPtr EC, TimePtr t,
                 const std::string &p)
  : m_impl(new Impl(c, sys, conf, EC, t, p)) {
  // empty
}

MPI_Comm Context::com() const {
  return m_impl->com;
}

Context::UnitsSystemPtr Context::unit_system() const {
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

const std::string& Context::prefix() const {
  return m_impl->prefix;
}

const Profiling& Context::profiling() const {
  return m_impl->profiling;
}

Context::Ptr context_from_options(MPI_Comm com, const std::string &prefix) {
  // unit system
  units::System::Ptr sys(new units::System);

  // configuration parameters
  Config::Ptr config = config_from_options(com, sys);
  print_config(3, com, *config);

  // time manager
  Time::Ptr time = time_from_options(com, config, sys);

  // enthalpy converter
  EnthalpyConverter::Ptr EC = enthalpy_converter_from_options(*config);

  return Context::Ptr(new Context(com, sys, config, EC, time, prefix));
}


} // end of namespace pism
