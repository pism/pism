/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#include <mpi.h>
#include <cmath>

#include "pism/util/io/File.hh"
#include "ConfigInterface.hh"
#include "Units.hh"
#include "pism_utilities.hh"
#include "pism_options.hh"
#include "error_handling.hh"

// include an implementation header so that we can allocate a DefaultConfig instance in
// config_from_options()
#include "Config.hh"
#include "pism/util/Logger.hh"

namespace pism {

struct Config::Impl {
  Impl(units::System::Ptr sys)
    : unit_system(sys) {
    // empty
  }

  units::System::Ptr unit_system;

  std::string filename;

  //! @brief Set of parameters set by the user. Used to warn about parameters that were set but were
  //! not used.
  std::set<std::string> parameters_set_by_user;
  //! @brief Set of parameters used in a run. Used to warn about parameters that were set but were
  //! not used.
  std::set<std::string> parameters_used;
};

Config::Config(units::System::Ptr system)
  : m_impl(new Impl(system)) {
  // empty
}

Config::~Config() {
  delete m_impl;
}

void Config::read(MPI_Comm com, const std::string &filename) {

  File file(com, filename, PISM_NETCDF3, PISM_READONLY); // OK to use netcdf3
  this->read(file);
}

void Config::read(const File &file) {
  this->read_impl(file);

  m_impl->filename = file.filename();
}

void Config::write(const File &file) const {
  this->write_impl(file);
}

void Config::write(MPI_Comm com, const std::string &filename, bool append) const {

  IO_Mode mode = append ? PISM_READWRITE : PISM_READWRITE_MOVE;

  File file(com, filename, PISM_NETCDF3, mode); // OK to use netcdf3

  this->write(file);
}

//! \brief Returns the name of the file used to initialize the database.
std::string Config::filename() const {
  return m_impl->filename;
}

void Config::import_from(const Config &other) {
  auto parameters = this->keys();

  for (auto p : other.all_doubles()) {
    if (member(p.first, parameters)) {
      this->set_numbers(p.first, p.second, CONFIG_USER);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unrecognized parameter %s in %s",
                                    p.first.c_str(), other.filename().c_str());
    }
  }

  for (auto p : other.all_strings()) {
    if (member(p.first, parameters)) {
      this->set_string(p.first, p.second, CONFIG_USER);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unrecognized parameter %s in %s",
                                    p.first.c_str(), other.filename().c_str());
    }
  }

  for (auto p : other.all_flags()) {
    if (member(p.first, parameters)) {
      this->set_flag(p.first, p.second, CONFIG_USER);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unrecognized parameter %s in %s",
                                    p.first.c_str(), other.filename().c_str());
    }
  }
}

const std::set<std::string>& Config::parameters_set_by_user() const {
  return m_impl->parameters_set_by_user;
}

const std::set<std::string>& Config::parameters_used() const {
  return m_impl->parameters_used;
}

bool Config::is_set(const std::string &name) const {
  return this->is_set_impl(name);
}

Config::Doubles Config::all_doubles() const {
  return this->all_doubles_impl();
}

double Config::get_number(const std::string &name, UseFlag flag) const {
  if (flag == REMEMBER_THIS_USE) {
    m_impl->parameters_used.insert(name);
  }
  return this->get_number_impl(name);
}

double Config::get_number(const std::string &name,
                          const std::string &units,
                          UseFlag flag) const {
  double value = this->get_number(name, flag);
  std::string input_units = this->units(name);

  try {
    return units::convert(m_impl->unit_system, value, input_units, units);
  } catch (RuntimeError &e) {
    e.add_context("converting \"%s\" from \"%s\" to \"%s\"",
                  name.c_str(), input_units.c_str(), units.c_str());
    throw;
  }
}

std::vector<double> Config::get_numbers(const std::string &name, UseFlag flag) const {
  if (flag == REMEMBER_THIS_USE) {
    m_impl->parameters_used.insert(name);
  }
  return this->get_numbers_impl(name);
}

std::vector<double> Config::get_numbers(const std::string &name,
                                        const std::string &units,
                                        UseFlag flag) const {
  auto value       = this->get_numbers(name, flag);
  auto input_units = this->units(name);

  try {
    units::Converter converter(m_impl->unit_system, input_units, units);
    for (unsigned int k = 0; k < value.size(); ++k) {
      value[k] = converter(value[k]);
    }
    return value;
  } catch (RuntimeError &e) {
    e.add_context("converting \"%s\" from \"%s\" to \"%s\"",
                  name.c_str(), input_units.c_str(), units.c_str());
    throw;
  }
}

void Config::set_number(const std::string &name, double value,
                        ConfigSettingFlag flag) {
  std::set<std::string> &set_by_user = m_impl->parameters_set_by_user;

  if (flag == CONFIG_USER) {
    set_by_user.insert(name);
  }

  // stop if we're setting the default value and this parameter was set by user already
  if (flag == CONFIG_DEFAULT and
      set_by_user.find(name) != set_by_user.end()) {
    return;
  }

  this->set_number_impl(name, value);
}

void Config::set_numbers(const std::string &name,
                         const std::vector<double> &values,
                         ConfigSettingFlag flag) {
  std::set<std::string> &set_by_user = m_impl->parameters_set_by_user;

  if (flag == CONFIG_USER) {
    set_by_user.insert(name);
  }

  // stop if we're setting the default value and this parameter was set by user already
  if (flag == CONFIG_DEFAULT and
      set_by_user.find(name) != set_by_user.end()) {
    return;
  }

  this->set_numbers_impl(name, values);
}

Config::Strings Config::all_strings() const {
  return this->all_strings_impl();
}

std::string Config::get_string(const std::string &name, UseFlag flag) const {
  if (flag == REMEMBER_THIS_USE) {
    m_impl->parameters_used.insert(name);
  }
  return this->get_string_impl(name);
}

void Config::set_string(const std::string &name,
                        const std::string &value,
                        ConfigSettingFlag flag) {
  std::set<std::string> &set_by_user = m_impl->parameters_set_by_user;

  if (flag == CONFIG_USER) {
    set_by_user.insert(name);
  }

  // stop if we're setting the default value and this parameter was set by user already
  if (flag == CONFIG_DEFAULT and
      set_by_user.find(name) != set_by_user.end()) {
    return;
  }

  this->set_string_impl(name, value);
}

Config::Flags Config::all_flags() const {
  return this->all_flags_impl();
}

bool Config::get_flag(const std::string& name, UseFlag flag) const {
  if (flag == REMEMBER_THIS_USE) {
    m_impl->parameters_used.insert(name);
  }
  return this->get_flag_impl(name);
}

void Config::set_flag(const std::string& name, bool value,
                         ConfigSettingFlag flag) {
  std::set<std::string> &set_by_user = m_impl->parameters_set_by_user;

  if (flag == CONFIG_USER) {
    set_by_user.insert(name);
  }

  // stop if we're setting the default value and this parameter was set by user already
  if (flag == CONFIG_DEFAULT and
      set_by_user.find(name) != set_by_user.end()) {
    return;
  }

  this->set_flag_impl(name, value);
}

static bool special_parameter(const std::string &name) {
  for (auto suffix : {"_doc", "_units", "_type", "_option", "_choices"}) {
    if (ends_with(name, suffix)) {
      return true;
    }
  }

  // The NetCDF-based configuration database stores parameters as attributes of a variable
  // and CF conventions require that all variables have a "long name."
  if (name == "long_name") {
    return true;
  }

  return false;
}

void print_config(const Logger &log, int verbosity_threshhold, const Config &config) {
  const int v = verbosity_threshhold;

  log.message(v,
             "### Strings:\n"
             "###\n");

  Config::Strings strings = config.all_strings();

  // find max. name size
  size_t max_name_size = 0;
  for (auto s : strings) {
    if (special_parameter(s.first)) {
      continue;
    }
    max_name_size = std::max(max_name_size, s.first.size());
  }

  // print strings
  for (auto s : strings) {
    std::string name  = s.first;
    std::string value = s.second;

    if (value.empty() or special_parameter(name)) {
      continue;
    }

    std::string padding(max_name_size - name.size(), ' ');

    if (config.type(name) == "keyword") {
      log.message(v, "  %s%s = \"%s\" (allowed choices: %s)\n",
                  name.c_str(), padding.c_str(), value.c_str(),
                  config.choices(name).c_str());
    } else {
      log.message(v, "  %s%s = \"%s\"\n", name.c_str(), padding.c_str(), value.c_str());
    }
  }

  log.message(v,
             "### Doubles:\n"
             "###\n");

  // find max. name size
  max_name_size = 0;
  for (auto d : config.all_doubles()) {
    max_name_size = std::max(max_name_size, d.first.size());
  }
  // print doubles
  for (auto d : config.all_doubles()) {
    std::string name  = d.first;
    double      value = d.second[0];

    std::string units = config.units(name); // will be empty if not set
    std::string padding(max_name_size - name.size(), ' ');

    if (fabs(value) >= 1.0e7 or fabs(value) <= 1.0e-4) {
      // use scientific notation if a number is big or small
      log.message(v, "  %s%s = %13.3e (%s)\n", name.c_str(), padding.c_str(), value, units.c_str());
    } else {
      log.message(v, "  %s%s = %13.5f (%s)\n", name.c_str(), padding.c_str(), value, units.c_str());
    }
  }

  log.message(v,
             "### Flags:\n"
             "###\n");

  // find max. name size
  max_name_size = 0;
  for (auto b : config.all_flags()) {
    max_name_size = std::max(max_name_size, b.first.size());
  }

  // print flags
  for (auto b : config.all_flags()) {
    std::string name  = b.first;
    std::string value = b.second ? "true" : "false";
    std::string padding(max_name_size - name.size(), ' ');

    log.message(v, "  %s%s = %s\n", name.c_str(), padding.c_str(), value.c_str());
  }

  log.message(v,
             "### List of configuration parameters ends here.\n"
             "###\n");
}

void print_unused_parameters(const Logger &log, int verbosity_threshhold,
                             const Config &config) {
  std::set<std::string> parameters_set = config.parameters_set_by_user();
  std::set<std::string> parameters_used = config.parameters_used();

  if (options::Bool("-options_left", "report unused options")) {
    verbosity_threshhold = log.get_threshold();
  }

  for (auto p : parameters_set) {

    if (special_parameter(p)) {
      continue;
    }

    if (parameters_used.find(p) == parameters_used.end()) {
      log.message(verbosity_threshhold,
                  "PISM WARNING: flag or parameter \"%s\" was set but was not used!\n",
                  p.c_str());

    }
  }
}

// command-line options

//! Get a flag from a command-line option.
/*!
 * Use the command-line option `option` to set the configuration parameter `parameter_name`.
 *
 * When called as `set_flag_from_option(config, "foo", "bar")`,
 *
 * sets the configuration parameter `bar` to `true` if
 *
 * - `-foo` is set (no argument)
 * - `-foo true` is set
 * - `-foo True` is set
 * - `-foo yes` is set
 *
 * sets the `bar` to `false` if
 *
 * - `-foo false` is set
 * - `-foo False` is set
 * - `-foo no` is set
 * - `-no_foo is set.
 *
 * `-foo X` with `X` not equal to `yes`, `no`, `true`, `True`, `false`, `False` results in an error.
 */
void set_flag_from_option(Config &config, const std::string &option,
                             const std::string &parameter_name) {

  // get the default value
  bool value = config.get_flag(parameter_name, Config::FORGET_THIS_USE);
  std::string doc = config.doc(parameter_name);

  // process the command-line option
  options::String opt("-" + option, doc, value ? "true" : "false", options::ALLOW_EMPTY);

  if (opt.is_set()) {
    if (opt.value() == ""     or
        opt.value() == "on"   or
        opt.value() == "yes"  or
        opt.value() == "true" or
        opt.value() == "True") { // Python's "True"

      value = true;

    } else if (opt.value() == "off"   or
               opt.value() == "no"    or
               opt.value() == "false" or
               opt.value() == "False") { // Python's "False"

      value = false;

    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s argument: %s",
                                    option.c_str(), opt.value().c_str());
    }
  }

  // For backward compatibility we allow disabling an option -foo by setting -no_foo.
  {
    bool no_foo_is_set = options::Bool("-no_" + option, doc);

    if (no_foo_is_set) {
      if (opt.is_set()) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Inconsistent command-line options:"
                                      " both -%s and -no_%s are set.\n",
                                      option.c_str(), option.c_str());
      } else {
        value = false;
      }
    }
  }

  config.set_flag(parameter_name, value, CONFIG_USER);
}

//! Sets a configuration parameter from a command-line option.
/*!
  If called as number_from_option("foo", "foo"), checks -foo and calls set("foo", value).

  Does nothing if -foo was not set.

  Note that no unit conversion is performed; parameters should be stored in
  input units and converted as needed. (This allows saving parameters without
  converting again.)
*/
void set_number_from_option(units::System::Ptr unit_system, Config &config, const std::string &name, const std::string &parameter) {
  options::Real option(unit_system,
                       "-" + name,
                       config.doc(parameter),
                       config.units(parameter),
                       config.get_number(parameter, Config::FORGET_THIS_USE));
  if (option.is_set()) {
    config.set_number(parameter, option, CONFIG_USER);
  }
}

/*!
 * Use a command-line option -option to set a parameter that is a list of numbers.
 *
 * The length of the list given as an argument to the command-line option has to be the
 * same as the length of the default value of the parameter *unless* the length of the
 * default value is less than 2. This default value length is used to disable this check.
 */
void set_number_list_from_option(Config &config, const std::string &option,
                                 const std::string &parameter) {
  auto default_value = config.get_numbers(parameter, Config::FORGET_THIS_USE);
  options::RealList list("-" + option, config.doc(parameter), default_value);

  if (list.is_set()) {
    if (default_value.size() < 2 and list->size() != default_value.size()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "Option -%s requires a list of %d numbers (got %d instead).",
                                    option.c_str(),
                                    (int)default_value.size(),
                                    (int)list->size());
    }

    config.set_numbers(parameter, list, CONFIG_USER);
  }
}

/*!
 * Use a command-line option -option to set a parameter that is a list of integers.
 *
 * The length of the list given as an argument to the command-line option has to be the
 * same as the length of the default value of the parameter *unless* the length of the
 * default value is less than 2. This default value length is used to disable this check.
 */
void set_integer_list_from_option(Config &config, const std::string &option,
                                  const std::string &parameter) {
  std::vector<int> default_value;

  for (auto v : config.get_numbers(parameter, Config::FORGET_THIS_USE)) {
    default_value.push_back(v);
  }

  options::IntegerList list("-" + option, config.doc(parameter), default_value);

  std::vector<double> value;
  for (auto v : list.value()) {
    value.push_back(v);
  }

  if (list.is_set()) {
    if (default_value.size() < 2 and value.size() != default_value.size()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "Option -%s requires a list of %d integers (got %d instead).",
                                    option.c_str(),
                                    (int)default_value.size(),
                                    (int)value.size());
    }

    config.set_numbers(parameter, value, CONFIG_USER);
  }
}

void set_integer_from_option(Config &config, const std::string &name, const std::string &parameter) {
  options::Integer option("-" + name, config.doc(parameter),
                          config.get_number(parameter, Config::FORGET_THIS_USE));
  if (option.is_set()) {
    config.set_number(parameter, option, CONFIG_USER);
  }
}

void set_string_from_option(Config &config, const std::string &name, const std::string &parameter) {

  options::String value("-" + name, config.doc(parameter),
                        config.get_string(parameter, Config::FORGET_THIS_USE));
  if (value.is_set()) {
    config.set_string(parameter, value, CONFIG_USER);
  }
}

//! \brief Set a keyword parameter from a command-line option.
/*!
 * This sets the parameter "parameter" after checking the "-name" command-line
 * option. This option requires an argument, which has to match one of the
 * keyword given in a comma-separated list "choices_list".
 */
void set_keyword_from_option(Config &config, const std::string &name,
                             const std::string &parameter,
                             const std::string &choices) {

  options::Keyword keyword("-" + name, config.doc(parameter), choices,
                           config.get_string(parameter, Config::FORGET_THIS_USE));

  if (keyword.is_set()) {
    config.set_string(parameter, keyword, CONFIG_USER);
  }
}

void set_parameter_from_options(units::System::Ptr unit_system,
                                Config &config, const std::string &name) {

  // skip special parameters ("attributes" of parameters)
  if (special_parameter(name)) {
    return;
  }

  // Use parameter name as its own command-line option by default. parameter_name_option can specify
  // a different (possibly shorter) command-line option.
  std::string option = name;

  if (not config.option(name).empty()) { // there is a short version of the command-line option
    std::string
      short_option = config.option(name),
      description  = config.doc(name);

    if (options::Bool("-" + short_option, description) or
        options::Bool("-no_" + short_option, description)) { // short option is set
      if (options::Bool("-" + option, description)) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "both -%s and -%s are set (please use one or the other)",
                                      option.c_str(), short_option.c_str());
      }

      // Use the short option only if the user set it, otherwise used the full (long) option below.
      option = short_option;
    }
  }

  std::string type = config.type(name);

  if (type == "string") {
    set_string_from_option(config, option, name);
  } else if (type == "flag") {
    set_flag_from_option(config, option, name);
  } else if (type == "number") {
    set_number_from_option(unit_system, config, option, name);
  } else if (type == "integer") {
    set_integer_from_option(config, option, name);
  } else if (type == "keyword") {
    set_keyword_from_option(config, option, name, config.choices(name));
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "parameter type \"%s\" is invalid", type.c_str());
  }
}

void set_config_from_options(units::System::Ptr unit_system,
                             Config &config) {
  for (const auto &d : config.all_doubles()) {
    set_parameter_from_options(unit_system, config, d.first);
  }

  for (const auto &s : config.all_strings()) {
    set_parameter_from_options(unit_system, config, s.first);
  }

  for (const auto &b : config.all_flags()) {
    set_parameter_from_options(unit_system, config, b.first);
  }

  // Energy modeling
  {
    options::Keyword energy("-energy",
                            "choose the energy model (one of 'none', 'cold', 'enthalpy')",
                            "none,cold,enthalpy", "enthalpy");

    if (energy.is_set()) {
      if (energy == "none") {
        config.set_flag("energy.enabled", false, CONFIG_USER);
        // Allow selecting cold ice flow laws in isothermal mode.
        config.set_flag("energy.temperature_based", true, CONFIG_USER);
      } else if (energy == "cold") {
        config.set_flag("energy.enabled", true, CONFIG_USER);
        config.set_flag("energy.temperature_based", true, CONFIG_USER);
      } else if (energy == "enthalpy") {
        config.set_flag("energy.enabled", true, CONFIG_USER);
        config.set_flag("energy.temperature_based", false, CONFIG_USER);
      } else {
        throw RuntimeError(PISM_ERROR_LOCATION, "this can't happen: options::Keyword validates input");
      }
    }
  }

  // -iterative_phi
   // 2,5,70,1,250,500,-300,700,1e-3

  {
    std::vector<double> defaults = {

      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_min"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_minup"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_max"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.dphi"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.dt"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.h_inv"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.topg_min"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.topg_max"),
      config.get_number("basal_yield_stress.mohr_coulomb.iterative_phi.dh_conv")
    };

    options::RealList iterative_phi("-iterative_phi", "phi_min, phi_minup, phi_max, dphi, dt_phi_inv, h_inv, topg_min, topg_max, dhdt_conv",
                                  defaults);
    if (iterative_phi.is_set()) {
      if (iterative_phi->size() != 9) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "option -iterative_phi expected 9 numbers; got %d",
                                      (int)iterative_phi->size());
      }

      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_min", iterative_phi[0]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_minup", iterative_phi[1]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.phi_max", iterative_phi[2]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.dphi", iterative_phi[3]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.dt", iterative_phi[4]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.h_inv", iterative_phi[5]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.topg_min", iterative_phi[6]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.topg_max", iterative_phi[7]);
      config.set_number("basal_yield_stress.mohr_coulomb.iterative_phi.dh_conv", iterative_phi[8]);

    }
  }


  // -topg_to_phi
  {
    std::vector<double> defaults = {
      config.get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min"),
      config.get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max"),
      config.get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min"),
      config.get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max")
    };

    options::RealList topg_to_phi("-topg_to_phi", "phi_min, phi_max, topg_min, topg_max",
                                  defaults);
    if (topg_to_phi.is_set()) {
      if (topg_to_phi->size() != 4) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "option -topg_to_phi expected 4 numbers; got %d",
                                      (int)topg_to_phi->size());
      }
      config.set_flag("basal_yield_stress.mohr_coulomb.topg_to_phi.enabled", true);
      config.set_number("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min", topg_to_phi[0]);
      config.set_number("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max", topg_to_phi[1]);
      config.set_number("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min", topg_to_phi[2]);
      config.set_number("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max", topg_to_phi[3]);
    }
  }
  // Ice shelves

  bool nu_bedrock = options::Bool("-nu_bedrock", "constant viscosity near margins");
  if (nu_bedrock) {
    config.set_flag("stress_balance.ssa.fd.lateral_drag.enabled", true, CONFIG_USER);
  }

  // Shortcuts

  // option "-pik" turns on a suite of PISMPIK effects (but NOT a calving choice,
  // and in particular NOT  "-calving eigen_calving")
  bool pik = options::Bool("-pik", "enable suite of PISM-PIK mechanisms");
  if (pik) {
    config.set_flag("stress_balance.calving_front_stress_bc", true, CONFIG_USER);
    config.set_flag("geometry.part_grid.enabled", true, CONFIG_USER);
    config.set_flag("geometry.remove_icebergs", true, CONFIG_USER);
    config.set_flag("geometry.grounded_cell_fraction", true, CONFIG_USER);
  }

  if (config.get_string("calving.methods").find("eigen_calving") != std::string::npos) {
    config.set_flag("geometry.part_grid.enabled", true, CONFIG_USER);
    // eigen-calving requires a wider stencil:
    config.set_number("grid.max_stencil_width", 3);
  }

  // all calving mechanisms require iceberg removal
  if (not config.get_string("calving.methods").empty()) {
    config.set_flag("geometry.remove_icebergs", true, CONFIG_USER);
  }

  // geometry.remove_icebergs requires part_grid
  if (config.get_flag("geometry.remove_icebergs")) {
    config.set_flag("geometry.part_grid.enabled", true, CONFIG_USER);
  }

  bool test_climate_models = options::Bool("-test_climate_models",
                                           "Disable ice dynamics to test climate models");
  if (test_climate_models) {
    config.set_string("stress_balance.model", "none", CONFIG_USER);
    config.set_flag("energy.enabled", false, CONFIG_USER);
    config.set_flag("age.enabled", false, CONFIG_USER);
    // let the user decide if they want to use "-no_mass" or not
  }

  // If frontal melt code includes floating ice, routing hydrology should include it also.
  if (config.get_string("hydrology.model") == "routing") {
    if (config.get_flag("frontal_melt.include_floating_ice")) {
      config.set_flag("hydrology.routing.include_floating_ice", true);
    }
  }

  if (config.get_flag("output.ISMIP6")) {
    // use MKS units in ISMIP6 mode
    config.set_flag("output.use_MKS", true);
  }
}

//! Create a configuration database using command-line options.
Config::Ptr config_from_options(MPI_Comm com, const Logger &log, units::System::Ptr sys) {

  DefaultConfig::Ptr config(new DefaultConfig(com, "pism_config", "-config", sys)),
    overrides(new DefaultConfig(com, "pism_overrides", "-config_override", sys));
  overrides->init(log);
  config->init_with_default(log);
  config->import_from(*overrides);
  set_config_from_options(sys, *config);

  return config;
}

ConfigWithPrefix::ConfigWithPrefix(Config::ConstPtr c, const std::string &prefix)
  : m_prefix(prefix), m_config(c) {
  // empty
}

double ConfigWithPrefix::get_number(const std::string &name) const {
  return m_config->get_number(m_prefix + name);
}

double ConfigWithPrefix::get_number(const std::string &name, const std::string &units) const {
  return m_config->get_number(m_prefix + name, units);
}

std::string ConfigWithPrefix::get_string(const std::string &name) const {
  return m_config->get_string(m_prefix + name);
}

bool ConfigWithPrefix::get_flag(const std::string& name) const {
  return m_config->get_flag(m_prefix + name);
}

void ConfigWithPrefix::reset_prefix(const std::string &prefix) {
  m_prefix = prefix;
}

std::set<std::string> Config::keys() const {
  std::set<std::string> result;

  for (auto p : all_doubles()) {
    result.insert(p.first);
  }

  for (auto p : all_strings()) {
    result.insert(p.first);
  }

  for (auto p : all_flags()) {
    result.insert(p.first);
  }

  return result;
}

std::string Config::doc(const std::string &parameter) const {
  return this->get_string(parameter + "_doc", Config::FORGET_THIS_USE);
}

std::string Config::units(const std::string &parameter) const {
  return this->get_string(parameter + "_units", Config::FORGET_THIS_USE);
}

std::string Config::type(const std::string &parameter) const {
  return this->get_string(parameter + "_type", Config::FORGET_THIS_USE);
}

std::string Config::option(const std::string &parameter) const {
  if (this->is_set(parameter + "_option")) {
    return this->get_string(parameter + "_option", Config::FORGET_THIS_USE);
  } else {
    return "";
  }
}

std::string Config::choices(const std::string &parameter) const {
  return this->get_string(parameter + "_choices", Config::FORGET_THIS_USE);
}



} // end of namespace pism
