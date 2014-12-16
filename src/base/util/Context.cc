/* Copyright (C) 2014 PISM Authors
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

Context::Context(Config &config)
  : m_config(config) {
}

MPI_Comm Context::com() const {
  return m_config.com();
}

Vars& Context::variables() {
  return m_variables;
}

const Vars& Context::variables() const {
  return m_variables;
}

Config &Context::config() {
  return m_config;
}


const Config& Context::config() const {
  return m_config;
}

double Context::convert(double value, const std::string &unit1, const std::string &unit2) const {
  return unit_system().convert(value, unit1, unit2);
}

UnitSystem Context::unit_system() const {
  return m_config.get_unit_system();
}

Time::Ptr Context::time() {
  return m_time;
}

const Time::Ptr Context::time() const {
  return m_time;
}

Profiling Context::profiling() const {
  return m_profiling;
}

} // end of namespace pism
