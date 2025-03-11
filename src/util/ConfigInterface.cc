/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025 PISM Authors
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
#include <cmath>                // std::round()
#include <cstdlib>              // realpath()


#include "pism/pism_config.hh"
#include "pism/util/io/File.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Units.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/IO_Flags.hh"

// include an implementation header so that we can allocate a NetCDFConfig instance in
// config_from_options()
#include "pism/util/NetCDFConfig.hh"
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

  File file(com, filename, io::PISM_NETCDF3, io::PISM_READONLY); // OK to use netcdf3
  this->read(file);
}

void Config::read(const File &file) {
  this->read_impl(file);

  m_impl->filename = file.name();
}

void Config::write(const File &file) const {
  this->write_impl(file);
}

void Config::write(MPI_Comm com, const std::string &filename, bool append) const {

  io::Mode mode = append ? io::PISM_READWRITE : io::PISM_READWRITE_MOVE;

  File file(com, filename, io::PISM_NETCDF3, mode); // OK to use netcdf3

  this->write(file);
}

//! \brief Returns the name of the file used to initialize the database.
std::string Config::filename() const {
  return m_impl->filename;
}

void Config::import_from(const Config &other) {
  auto parameters = this->keys();

  for (const auto &p : other.all_doubles()) {
    if (member(p.first, parameters)) {
      this->set_numbers(p.first, p.second, CONFIG_USER);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unrecognized parameter %s in %s",
                                    p.first.c_str(), other.filename().c_str());
    }
  }

  for (const auto &p : other.all_strings()) {
    if (member(p.first, parameters)) {
      this->set_string(p.first, p.second, CONFIG_USER);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unrecognized parameter %s in %s",
                                    p.first.c_str(), other.filename().c_str());
    }
  }

  for (const auto &p : other.all_flags()) {
    if (member(p.first, parameters)) {
      this->set_flag(p.first, p.second, CONFIG_USER);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unrecognized parameter %s in %s",
                                    p.first.c_str(), other.filename().c_str());
    }
  }
}

void Config::resolve_filenames() {
  for (const auto &s : all_strings()) {
    auto parameter = s.first;
    auto value = s.second;

    if (value.empty()) {
      continue;
    }

    auto last_token = pism::split(parameter, '.').back();

    if (last_token == "file") {
      char *resolved_path = realpath(value.c_str(), NULL);

      if (resolved_path != NULL) {
        set_string(parameter, resolved_path, CONFIG_USER);
        free(resolved_path);
      }
      // Note: we keep the old value if `realpath()` failed
    }
  } // end of the loop over all strings
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

bool Config::is_valid_number(const std::string &name) const {
  auto value = get_number(name, FORGET_THIS_USE);
  auto min = valid_min(name);

  if (std::get<0>(min)) {
    auto valid_min = std::get<1>(min);

    if (value < valid_min) {
      return false;
    }
  }

  auto max = valid_max(name);
  if (std::get<0>(max)) {
    auto valid_max = std::get<1>(max);

    if (value > valid_max) {
      return false;
    }
  }

  return true;
}

double Config::get_number(const std::string &name, UseFlag flag) const {
  auto value = get_number_impl(name);

  if (flag == REMEMBER_THIS_USE) {
    // check the valid range (if set) and remember that this parameter was used
    //
    // note that we don't check the valid range when flag == FORGET_THIS_USE. This way we
    // can get the default value of a parameter. Parameters without a default value should
    // be set to values outside of their respective valid ranges (if possible).
    m_impl->parameters_used.insert(name);

    if (type(name) == "integer" and std::round(value) != value) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "integer parameter '%s' was set to a number with a non-zero fractional part (%f)",
          name.c_str(), value);
    }

    auto min = valid_min(name);
    if (std::get<0>(min)) {
      auto valid_min = std::get<1>(min);

      if (value < valid_min) {
        if (type(name) == "integer") {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "Please set '%s' to a number greater than or equal to %d",
                                        name.c_str(), (int)valid_min);
        }
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Please set '%s' to a number greater than or equal to %f",
                                      name.c_str(), valid_min);
      }
    }

    auto max = valid_max(name);
    if (std::get<0>(max)) {
      auto valid_max = std::get<1>(max);

      if (value > valid_max) {
        if (type(name) == "integer") {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "Please set '%s' to a number less than or equal to %d",
                                        name.c_str(), (int)valid_max);
        }
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Please set '%s' to a number less than or equal to %f",
                                      name.c_str(), valid_max);
      }
    }
  }

  return value;
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
  for (const auto &suffix : {"_doc", "_units", "_type", "_option", "_choices", "_valid_min", "_valid_max"}) {
    if (ends_with(name, suffix)) {
      return true;
    }
  }

  // The NetCDF-based configuration database stores parameters as attributes of a variable
  // and CF conventions require that all variables have a "long name."
  return (name == "long_name");
}

void print_config(const Logger &log, int verbosity_threshhold, const Config &config) {
  const int v = verbosity_threshhold;

  log.message(v,
             "### Strings:\n"
             "###\n");

  Config::Strings strings = config.all_strings();

  // find max. name size
  size_t max_name_size = 0;
  for (const auto &s : strings) {
    if (special_parameter(s.first)) {
      continue;
    }
    max_name_size = std::max(max_name_size, s.first.size());
  }

  // print strings
  for (const auto &s : strings) {
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
  for (const auto &d : config.all_doubles()) {
    max_name_size = std::max(max_name_size, d.first.size());
  }
  // print doubles
  for (auto d : config.all_doubles()) {
    std::string name  = d.first;
    double      value = d.second[0];

    if (special_parameter(name)) {
      continue;
    }

    std::string units = config.units(name); // will be empty if not set
    std::string padding(max_name_size - name.size(), ' ');

    const double large = 1.0e7;
    const double small = 1.0e-4;
    if (fabs(value) >= large or fabs(value) <= small) {
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
  for (const auto &b : config.all_flags()) {
    max_name_size = std::max(max_name_size, b.first.size());
  }

  // print flags
  for (const auto &b : config.all_flags()) {
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

  for (const auto &p : parameters_set) {

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
  bool value      = config.get_flag(parameter_name, Config::FORGET_THIS_USE);
  std::string doc = config.doc(parameter_name);

  // process the command-line option
  options::String opt("-" + option, doc, value ? "true" : "false", options::ALLOW_EMPTY);

  if (opt.is_set()) {
    if (member(opt.value(), { "", "on", "yes", "true", "True" })) {

      value = true;

    } else if (member(opt.value(), { "off", "no", "false", "False" })) {

      value = false;

    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s argument: %s", option.c_str(),
                                    opt.value().c_str());
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
      }

      value = false;
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
void set_number_from_option(units::System::Ptr unit_system, Config &config, const std::string &option,
                            const std::string &parameter) {
  options::Real opt(unit_system, "-" + option, config.doc(parameter), config.units(parameter),
                    config.get_number(parameter, Config::FORGET_THIS_USE));
  if (opt.is_set()) {
    config.set_number(parameter, opt, CONFIG_USER);
  }
}

void set_integer_from_option(Config &config, const std::string &option,
                             const std::string &parameter) {
  options::Integer opt("-" + option, config.doc(parameter),
                       (int)config.get_number(parameter, Config::FORGET_THIS_USE));
  if (opt.is_set()) {
    config.set_number(parameter, opt, CONFIG_USER);
  }
}

void set_string_from_option(Config &config, const std::string &option, const std::string &parameter) {

  options::String value("-" + option, config.doc(parameter),
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
void set_keyword_from_option(Config &config, const std::string &option, const std::string &parameter,
                             const std::string &choices) {

  options::Keyword keyword("-" + option, config.doc(parameter), choices,
                           config.get_string(parameter, Config::FORGET_THIS_USE));

  if (keyword.is_set()) {
    config.set_string(parameter, keyword, CONFIG_USER);
  }
}

void set_parameter_from_options(units::System::Ptr unit_system, Config &config,
                                const std::string &name) {

  // skip special parameters ("attributes" of parameters)
  if (special_parameter(name)) {
    return;
  }

  // Use parameter name as its own command-line option by default. parameter_name_option can specify
  // a different (possibly shorter) command-line option.
  std::string option = name;

  if (not config.option(name).empty()) { // there is a short version of the command-line option
    std::string short_option = config.option(name);
    std::string description  = config.doc(name);

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
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "parameter type \"%s\" is invalid",
                                  type.c_str());
  }
}

void set_config_from_options(units::System::Ptr unit_system, Config &config) {
  for (const auto &d : config.all_doubles()) {
    set_parameter_from_options(unit_system, config, d.first);
  }

  for (const auto &s : config.all_strings()) {
    set_parameter_from_options(unit_system, config, s.first);
  }

  for (const auto &b : config.all_flags()) {
    set_parameter_from_options(unit_system, config, b.first);
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

  // Special command-line options for "-surface elevation,...":
  {
    options::String T("-ice_surface_temp", "ice surface temperature parameterization");

    if (T.is_set()) {
      auto IST = parse_number_list(T);

      if (IST.size() != 4) {
        throw RuntimeError(PISM_ERROR_LOCATION, "option -ice_surface_temp requires an argument"
                                                " (comma-separated list of 4 numbers)");
      }

      config.set_number("surface.elevation_dependent.T_min", IST[0]);
      config.set_number("surface.elevation_dependent.T_max", IST[1]);
      config.set_number("surface.elevation_dependent.z_T_min", IST[2]);
      config.set_number("surface.elevation_dependent.z_T_max", IST[3]);
    }

    options::String M("-climatic_mass_balance", "climatic mass balance parameterization");

    if (M.is_set()) {
      auto CMB = parse_number_list(M);

      if (CMB.size() != 5) {
        throw RuntimeError(PISM_ERROR_LOCATION, "-climatic_mass_balance requires an argument"
                                                " (comma-separated list of 5 numbers)");
      }

      config.set_number("surface.elevation_dependent.M_min", CMB[0]);
      config.set_number("surface.elevation_dependent.M_max", CMB[1]);
      config.set_number("surface.elevation_dependent.z_M_min", CMB[2]);
      config.set_number("surface.elevation_dependent.z_ELA", CMB[3]);
      config.set_number("surface.elevation_dependent.z_M_max", CMB[4]);
    }

    options::String limits("-climatic_mass_balance_limits",
                           "lower and upper limits of the climatic mass balance");

    if (limits.is_set()) {

      auto L = parse_number_list(limits);

      if (L.size() != 2) {
        throw RuntimeError(PISM_ERROR_LOCATION, "-climatic_mass_balance_limits requires an argument"
                                                " (a comma-separated list of 2 numbers)");
      }

      units::Converter meter_per_second(unit_system, "m year^-1", "m second^-1");

      config.set_number("surface.elevation_dependent.M_limit_min", meter_per_second(L[0]));
      config.set_number("surface.elevation_dependent.M_limit_max", meter_per_second(L[1]));
    }
  }
}

//! Create a configuration database using command-line options.
Config::Ptr config_from_options(MPI_Comm com, const Logger &log, units::System::Ptr unit_system) {

  using T         = NetCDFConfig;
  auto  config    = std::make_shared<T>(com, "pism_config", unit_system);
  auto  overrides = std::make_shared<T>(com, "pism_overrides", unit_system);

  options::String config_filename("-config", "Config file name", pism::config_file);
  options::String override_filename("-config_override", "Config override file name");

  // config_filename is always non-empty because it has a default value
  // (pism::config_file).
  config->read(com, config_filename);

  // override_filename may be empty if -config_override was not set
  if (override_filename.is_set()) {
    overrides->read(com, override_filename);
    config->import_from(*overrides);
  }

  set_config_from_options(unit_system, *config);

  config->resolve_filenames();

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

bool ConfigWithPrefix::get_flag(const std::string &name) const {
  return m_config->get_flag(m_prefix + name);
}

void ConfigWithPrefix::reset_prefix(const std::string &prefix) {
  m_prefix = prefix;
}

std::set<std::string> Config::keys() const {
  std::set<std::string> result;

  for (const auto &p : all_doubles()) {
    result.insert(p.first);
  }

  for (const auto &p : all_strings()) {
    result.insert(p.first);
  }

  for (const auto &p : all_flags()) {
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
  }

  return "";
}

std::string Config::choices(const std::string &parameter) const {
  return this->get_string(parameter + "_choices", Config::FORGET_THIS_USE);
}

std::pair<bool, double> Config::valid_min(const std::string &parameter) const {
  if (is_set(parameter + "_valid_min")) {
    return { true, get_number(parameter + "_valid_min", Config::FORGET_THIS_USE) };
  }
  return { false, {} };
}

std::pair<bool, double> Config::valid_max(const std::string &parameter) const {
  if (is_set(parameter + "_valid_max")) {
    return { true, get_number(parameter + "_valid_max", Config::FORGET_THIS_USE) };
  }
  return { false, {} };
}

} // end of namespace pism
