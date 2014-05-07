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

#include "PISMConfig.hh"
#include "pism_options.hh"
#include <sstream>

namespace pism {

PISMConfig::PISMConfig(MPI_Comm com, std::string name, PISMUnitSystem unit_system)
  : m_com(com),
    m_unit_system(unit_system),
    m_data(name, unit_system) {
  m_options_left_set = false;
  PISMOptionsIsSet("-options_left", m_options_left_set);
  m_unit_system = unit_system;
}

PISMConfig::~PISMConfig() {
  warn_about_unused_parameters();
}

void PISMConfig::set_string(const std::string &name, const std::string &value) {
  m_data.set_string(name, value);
}

void PISMConfig::set_double(const std::string &name, double value) {
  m_data.set_double(name, value);
}

bool PISMConfig::is_set(const std::string &name) const {
  return m_data.has_attribute(name);
}

PetscErrorCode PISMConfig::read(std::string filename) {
  PetscErrorCode ierr;

  PIO nc(m_com, "netcdf3", m_unit_system); // OK to use netcdf3

  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);

  ierr = this->read(nc); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMConfig::write(std::string filename, bool append) const {
  PetscErrorCode ierr;

  PIO nc(m_com, "netcdf3", m_unit_system); // OK to use netcdf3

  PISM_IO_Mode mode = PISM_READWRITE;
  if (append == false) {
    mode = PISM_READWRITE_MOVE;
  }

  ierr = nc.open(filename, mode); CHKERRQ(ierr);

  ierr = this->write(nc); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Read boolean flags and double parameters from a NetCDF file.
/*!
  Erases all the present parameters before reading.
*/
PetscErrorCode PISMConfig::read(const PIO &nc) {

  PetscErrorCode ierr = nc.read_attributes(m_data.get_name(), m_data); CHKERRQ(ierr);

  m_config_filename = nc.inq_filename();

  return 0;
}

//! Write a config variable to a file (with all its attributes).
PetscErrorCode PISMConfig::write(const PIO &nc) const {
  PetscErrorCode ierr;
  bool variable_exists;

  ierr = nc.inq_var(m_data.get_name(), variable_exists); CHKERRQ(ierr);

  if (variable_exists == false) {
    ierr = nc.redef(); CHKERRQ(ierr);

    ierr = nc.def_var(m_data.get_name(),
                      PISM_BYTE, std::vector<std::string>()); CHKERRQ(ierr);

    ierr = nc.write_attributes(m_data, PISM_DOUBLE, false); CHKERRQ(ierr);
  } else {
    ierr = nc.write_attributes(m_data, PISM_DOUBLE, false); CHKERRQ(ierr);
  }

  return 0;
}

double PISMConfig::get_quiet(std::string name) const {
  const NCVariable::DoubleAttrs& doubles = m_data.get_all_doubles();
  if (doubles.find(name) != doubles.end()) {
    return m_data.get_double(name);
  } else {
    PetscPrintf(m_com, "PISM ERROR: parameter '%s' is unset. (Parameters read from '%s'.)\n",
                name.c_str(), m_config_filename.c_str());
    PISMEnd();
  }

  return 0;                     // can't happen
}

std::string PISMConfig::get_string_quiet(std::string name) const {
  const NCVariable::StringAttrs& strings = m_data.get_all_strings();
  if (strings.find(name) != strings.end())
    return m_data.get_string(name);
  else {
    PetscPrintf(m_com,
                "PISM ERROR: Parameter '%s' was not set. (Read from '%s'.)\n",
                name.c_str(), m_config_filename.c_str());
    PISMEnd();
  }

  return std::string();         // will never happen
}

bool PISMConfig::get_flag_quiet(std::string name) const {
  const NCVariable::StringAttrs& strings = m_data.get_all_strings();
  NCVariable::StringAttrs::const_iterator j = strings.find(name);
  if (j != strings.end()) {

    const std::string &value = j->second;

    if ((value == "false") ||
        (value == "no") ||
        (value == "off"))
      return false;

    if ((value == "true") ||
        (value == "yes") ||
        (value == "on"))
      return true;

    PetscPrintf(m_com,
                "PISM ERROR: Parameter '%s' (%s) cannot be interpreted as a boolean.\n"
                "            Please make sure that it is set to one of 'true', 'yes', 'on', 'false', 'no', 'off'.\n",
                name.c_str(), value.c_str());
    PISMEnd();
  }

  PetscPrintf(m_com,
              "PISM ERROR: Parameter '%s' was not set. (Read from '%s'.)\n",
              name.c_str(), m_config_filename.c_str());
  PISMEnd();

  return true;                  // will never happen
}


//! Returns a `double` parameter. Stops if it was not found.
double PISMConfig::get(const std::string &name) const {
  if (m_options_left_set)
    m_parameters_used.insert(name);

  return this->get_quiet(name);
}

double PISMConfig::get(std::string name, std::string u1, std::string u2) const {
  // always use get() (*not* _quiet) here
  return m_unit_system.convert(this->get(name),  u1.c_str(),  u2.c_str());
}


//! Returns a boolean flag by name. Unset flags are treated as if they are set to 'false'.
/*!
  Strings "false", "no", "off" are interpreted as 'false'; "true", "on", "yes" -- as 'true'.

  Any other string produces an error.
*/
bool PISMConfig::get_flag(const std::string &name) const {
  if (m_options_left_set)
    m_parameters_used.insert(name);

  return this->get_flag_quiet(name);
}

//! \brief Get a string attribute by name.
std::string PISMConfig::get_string(const std::string &name) const {
  if (m_options_left_set)
    m_parameters_used.insert(name);

  return this->get_string_quiet(name);
}

//! Set a value of a boolean flag.
void PISMConfig::set_flag(const std::string &name, bool value) {
  if (value)
    m_data.set_string(name, "true");
  else
    m_data.set_string(name, "false");
}

//! Get a flag from a command-line option.
/*!
  If called as flag_from_option("foo", "foo"), checks both -foo and -no_foo.

  \li if -foo is set, calls set_flag("foo", true),

  \li if -no_foo is set, calls set_flag("foo", false),

  \li if both are set, prints an error message and stops,

  \li if none, does nothing.

*/
PetscErrorCode PISMConfig::flag_from_option(std::string name, std::string flag) {
  PetscErrorCode ierr;
  bool foo = false,
    no_foo = false;

  ierr = PISMOptionsIsSet("-" + name, get_string_quiet(flag + "_doc"), foo); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-no_" + name, no_foo); CHKERRQ(ierr);

  if (foo && no_foo) {
    PetscPrintf(m_com, "PISM ERROR: Inconsistent command-line options: both -%s and -no_%s are set.\n",
                name.c_str(), name.c_str());
    PISMEnd();
  }

  if (foo)
    set_flag_from_option(flag, true);

  if (no_foo)
    set_flag_from_option(flag, false);

  return 0;
}

//! Sets a configuration parameter from a command-line option.
/*!
  If called as scalar_from_option("foo", "foo"), checks -foo and calls set("foo", value).

  Does nothing if -foo was not set.

  Note that no unit conversion is performed; parameters should be stored in
  input units and converted as needed. (This allows saving parameters without
  converting again.)
*/
PetscErrorCode PISMConfig::scalar_from_option(std::string name, std::string parameter) {
  PetscErrorCode ierr;
  double value = get_quiet(parameter);
  bool flag;

  ierr = PISMOptionsReal("-" + name,
                         get_string_quiet(parameter + "_doc"),
                         value, flag); CHKERRQ(ierr);
  if (flag) {
    this->set_scalar_from_option(parameter, value);
  }

  return 0;
}

PetscErrorCode PISMConfig::string_from_option(std::string name, std::string parameter) {
  PetscErrorCode ierr;
  std::string value = get_string_quiet(parameter);
  bool flag;

  ierr = PISMOptionsString("-" + name,
                           get_string_quiet(parameter + "_doc"),
                           value, flag); CHKERRQ(ierr);
  if (flag) {
    this->set_string_from_option(parameter, value);
  }

  return 0;
}

//! \brief Set a keyword parameter from a command-line option.
/*!
 * This sets the parameter "parameter" after checking the "-name" command-line
 * option. This option requires an argument, which has to match one of the
 * keyword given in a comma-separated list "choices_list".
 */
PetscErrorCode PISMConfig::keyword_from_option(std::string name,
                                               std::string parameter,
                                               std::string choices_list) {
  PetscErrorCode ierr;
  std::istringstream arg(choices_list);
  std::set<std::string> choices;
  std::string keyword, tmp;
  bool flag;

  // Split the list:
  while (getline(arg, tmp, ','))
    choices.insert(tmp);

  ierr = PISMOptionsList(m_com, "-" + name,
                         get_string_quiet(parameter + "_doc"),
                         choices,
                         get_string_quiet(parameter), keyword, flag); CHKERRQ(ierr);

  if (flag) {
    this->set_string_from_option(parameter, keyword);
  }

  return 0;
}

PetscErrorCode PISMConfig::set_flag_from_option(std::string name, bool value) {

  m_parameters_set.insert(name);

  this->set_flag(name, value);

  return 0;
}

PetscErrorCode PISMConfig::set_scalar_from_option(std::string name, double value) {

  m_parameters_set.insert(name);

  m_data.set_double(name, value);

  return 0;
}

PetscErrorCode PISMConfig::set_string_from_option(std::string name, std::string value) {

  m_parameters_set.insert(name);

  m_data.set_string(name, value);

  return 0;
}

PetscErrorCode PISMConfig::set_keyword_from_option(std::string name, std::string value) {

  this->set_string_from_option(name, value);

  return 0;
}


//! Print all the attributes of a configuration variable.
PetscErrorCode PISMConfig::print_to_stdout(int vt) const {
  PetscErrorCode ierr;

  ierr = verbPrintf(vt, m_com, "PISM parameters read from %s:\n",
                    m_config_filename.c_str());

  ierr = m_data.report_to_stdout(m_com, vt); CHKERRQ(ierr);

  std::set<std::string>::const_iterator k;
  std::string output;

  for (k = m_parameters_set.begin(); k != m_parameters_set.end(); ++k) {

    if (ends_with(*k, "_doc"))
      continue;

    if (k == m_parameters_set.begin())
      output += *k;
    else
      output += std::string(", ") + (*k);
  }

  if (output.empty() == false) {
    ierr = verbPrintf(vt, m_com, "PISM flags and parameters set from the command line:\n  %s\n",
                      output.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Returns the name of the file used to initialize the database.
std::string PISMConfig::get_config_filename() const {
  return m_config_filename;
}

PISMUnitSystem PISMConfig::get_unit_system() const {
  return m_unit_system;
}

const NCVariable& PISMConfig::get_data() const {
  return m_data;
}

//! Imports values from the other config variable, silently overwriting present values.
void PISMConfig::import_from(const PISMConfig &other) {
  const NCVariable::DoubleAttrs &other_doubles = other.get_data().get_all_doubles();
  NCVariable::DoubleAttrs::const_iterator j;
  for (j = other_doubles.begin(); j != other_doubles.end(); ++j) {
    m_data.set_doubles(j->first, j->second);
    m_parameters_set.insert(j->first);
  }

  const NCVariable::StringAttrs &other_strings = other.get_data().get_all_strings();
  NCVariable::StringAttrs::const_iterator k;
  for (k = other_strings.begin(); k != other_strings.end(); ++k) {
    m_data.set_string(k->first, k->second);
    m_parameters_set.insert(k->first);
  }
}

//! Update values from the other config variable, overwriting present values but avoiding adding new ones.
void PISMConfig::update_from(const PISMConfig &other) {
  const NCVariable::DoubleAttrs &doubles = m_data.get_all_doubles();
  const NCVariable::DoubleAttrs &other_doubles = other.get_data().get_all_doubles();
  NCVariable::DoubleAttrs::const_iterator i, j;

  for (j = doubles.begin(); j != doubles.end(); ++j) {
    i = other_doubles.find(j->first);
    if (i != other_doubles.end()) {
      m_data.set_doubles(j->first, i->second);
    }
  }

  const NCVariable::StringAttrs &strings = m_data.get_all_strings();
  const NCVariable::StringAttrs &other_strings = other.get_data().get_all_strings();
  NCVariable::StringAttrs::const_iterator k, m;

  for (k = strings.begin(); k != strings.end(); ++k) {
    m = other_strings.find(k->first);
    if (m != other_strings.end()) {
      m_data.set_string(k->first, m->second);
    }
  }
}

PetscErrorCode PISMConfig::warn_about_unused_parameters() const {
  PetscErrorCode ierr;

  if (m_options_left_set == false)
    return 0;

  std::set<std::string>::const_iterator k;
  for (k = m_parameters_set.begin(); k != m_parameters_set.end(); ++k) {

    if (ends_with(*k, "_doc"))
      continue;

    if (m_parameters_used.find(*k) == m_parameters_used.end()) {
      ierr = verbPrintf(2, m_com,
                        "PISM WARNING: flag or parameter \"%s\" was set but was not used!\n",
                        k->c_str()); CHKERRQ(ierr);

    }
  }

  return 0;
}

} // end of namespace pism
