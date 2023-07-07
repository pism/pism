/* Copyright (C) 2014, 2015, 2016, 2017, 2019, 2021, 2023 PISM Authors
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

#include "pism/util/Config.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/pism_config.hh"  // pism::config_file
#include "pism/util/io/IO_Flags.hh"

namespace pism {

NetCDFConfig::NetCDFConfig(MPI_Comm com, const std::string &name, units::System::Ptr system)
  : Config(system),
    m_com(com),
    m_data(name, system) {
}

NetCDFConfig::~NetCDFConfig() {
}

bool NetCDFConfig::is_set_impl(const std::string &name) const {
  return m_data.has_attribute(name);
}

// doubles

double NetCDFConfig::get_number_impl(const std::string &name) const {
  const VariableMetadata::DoubleAttrs &doubles = m_data.all_doubles();
  if (doubles.find(name) != doubles.end()) {
    return m_data.get_number(name);
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "parameter '%s' is unset. (Parameters read from '%s'.)",
                                name.c_str(), m_config_filename.c_str());

  return 0; // can't happen
}

std::vector<double> NetCDFConfig::get_numbers_impl(const std::string &name) const {
  const VariableMetadata::DoubleAttrs& doubles = m_data.all_doubles();
  if (doubles.find(name) != doubles.end()) {
    return m_data.get_numbers(name);
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "parameter '%s' is unset. (Parameters read from '%s'.)",
                                name.c_str(), m_config_filename.c_str());

  return {};                    // can't happen
}

Config::Doubles NetCDFConfig::all_doubles_impl() const {
  Doubles result;

  for (const auto& d : m_data.all_doubles()) {
    result[d.first] = d.second;
  }
  return result;
}


void NetCDFConfig::set_number_impl(const std::string &name, double value) {
  m_data.set_number(name, value);
}

void NetCDFConfig::set_numbers_impl(const std::string &name,
                                    const std::vector<double> &values) {
  m_data.set_numbers(name, values);
}

// strings

std::string NetCDFConfig::get_string_impl(const std::string &name) const {
  const VariableMetadata::StringAttrs& strings = m_data.all_strings();
  if (strings.find(name) != strings.end()) {
    return m_data[name];
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "Parameter '%s' was not set. (Read from '%s'.)\n", name.c_str(),
                                m_config_filename.c_str());

  return std::string();         // will never happen
}

Config::Strings NetCDFConfig::all_strings_impl() const {
  VariableMetadata::StringAttrs strings = m_data.all_strings();
  Strings result;

  for (const auto& s : strings) {
    std::string name = s.first;
    std::string value = s.second;

    auto k = strings.find(name + "_type");
    if (k != strings.end() and k->second == "flag") {
      // Flags are stored as strings. Skip them.
      continue;
    }

    result[name] = value;
  }
  return result;
}

void NetCDFConfig::set_string_impl(const std::string &name, const std::string &value) {
  m_data[name] = value;
}

// flags

static bool string_is_false(const std::string &value) {
  return value == "false" or value == "off" or value == "no";
}

static bool string_is_true(const std::string &value) {
  return value == "true" or value == "on" or value == "yes";
}

bool NetCDFConfig::get_flag_impl(const std::string &name) const {
  const VariableMetadata::StringAttrs& strings = m_data.all_strings();
  auto j = strings.find(name);
  if (j != strings.end()) {

    const std::string &value = j->second;

    if (string_is_false(value)) {
      return false;
    }

    if (string_is_true(value)) {
      return true;
    }

    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Parameter '%s' (%s) cannot be interpreted as a flag.\n"
                                  "Please make sure that it is set to one of 'true', 'yes', 'on', 'false', 'no', 'off'.",
                                  name.c_str(), value.c_str());
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Parameter '%s' was not set. (Read from '%s'.)",
                                name.c_str(), m_config_filename.c_str());

  return true;                  // will never happen
}

Config::Flags NetCDFConfig::all_flags_impl() const {
  Flags result;

  for (const auto& b : m_data.all_strings()) {
    std::string name = b.first;
    std::string value = b.second;

    if (string_is_true(value)) {
      result[name] = true;
    } else if (string_is_false(value)) {
      result[name] = false;
    }
  }
  return result;
}

//! Set a value of a flag flag.
void NetCDFConfig::set_flag_impl(const std::string &name, bool value) {
  m_data[name] = value ? "true" : "false";
}

// file I/O

//! Read flag flags and double parameters from a NetCDF file.
/*!
  Erases all the present parameters before reading.
*/
void NetCDFConfig::read_impl(const File &nc) {

  io::read_attributes(nc, m_data.get_name(), m_data);

  m_config_filename = nc.filename();
}

//! Write a config variable to a file (with all its attributes).
void NetCDFConfig::write_impl(const File &nc) const {

  bool variable_exists = nc.find_variable(m_data.get_name());

  if (not variable_exists) {
    nc.define_variable(m_data.get_name(), io::PISM_BYTE, {});

    io::write_attributes(nc, m_data, io::PISM_DOUBLE);
  } else {
    io::write_attributes(nc, m_data, io::PISM_DOUBLE);
  }
}


//! Config that respects command-line options and stores data in a NetCDF variable.
DefaultConfig::DefaultConfig(MPI_Comm com,
                             const std::string &variable_name,
                             const std::string &option,
                             units::System::Ptr system)
  : NetCDFConfig(com, variable_name, system),
    m_option(option) {
  // empty
}

void DefaultConfig::init(const Logger &log, bool use_default_path) {
  options::String file(m_option,
                       "Name of the file to read " + m_data.get_name() + " from",
                       pism::config_file);
  if (use_default_path or file.is_set()) {
    this->read(m_com, file);
    log.message(2, "Reading configuration parameters (%s) from file '%s'.\n",
                m_data.get_name().c_str(), file->c_str());
  }
}

void DefaultConfig::init_with_default(const Logger &log) {
  this->init(log, true);
}

void DefaultConfig::init(const Logger &log) {
  this->init(log, false);
}

} // end of namespace pism
