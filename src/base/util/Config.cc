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

#include <string>
#include "Config.hh"

namespace pism {

Config::~Config() {
  // empty; required because this class has virtual methods
}

void Config::set_double(const std::string &name, double value) {
  this->set_double_impl(name, value);
}

void Config::set_boolean(const std::string &name, bool value) {
  this->set_boolean_impl(name, value);
}

void Config::set_string(const std::string &name, const std::string &value) {
  this->set_string_impl(name, value);
}

double Config::get_double(const std::string &name) const {
  return this->get_double_impl(name);
}

double Config::get(const std::string &name) const {
  return this->get_double(name);
}

bool Config::get_boolean(const std::string &name) const {
  return this->get_boolean_impl(name);
}

std::string Config::get_string(const std::string &name) const {
  return this->get_string_impl(name);
}

} // end of namespace pism
