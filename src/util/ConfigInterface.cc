/* Copyright (C) 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/util/io/PIO.hh"
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

void Config::read(MPI_Comm com, const std::string &file) {

  PIO nc(com, "netcdf3", file, PISM_READONLY); // OK to use netcdf3
  this->read(nc);
}

void Config::read(const PIO &nc) {
  this->read_impl(nc);

  m_impl->filename = nc.inq_filename();
}

void Config::write(const PIO &nc) const {
  this->write_impl(nc);
}

void Config::write(MPI_Comm com, const std::string &file, bool append) const {

  IO_Mode mode = append ? PISM_READWRITE : PISM_READWRITE_MOVE;

  PIO nc(com, "netcdf3", file, mode); // OK to use netcdf3

  this->write(nc);
}

//! \brief Returns the name of the file used to initialize the database.
std::string Config::filename() const {
  return m_impl->filename;
}

void Config::import_from(const Config &other) {
  auto parameters = this->keys();

  for (auto p : other.all_doubles()) {
    if (member(p.first, parameters)) {
      this->set_double(p.first, p.second, CONFIG_USER);
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

  for (auto p : other.all_booleans()) {
    if (member(p.first, parameters)) {
      this->set_boolean(p.first, p.second, CONFIG_USER);
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

double Config::get_double(const std::string &name, UseFlag flag) const {
  if (flag == REMEMBER_THIS_USE) {
    m_impl->parameters_used.insert(name);
  }
  return this->get_double_impl(name);
}

double Config::get_double(const std::string &name,
                          const std::string &units,
                          UseFlag flag) const {
  double value = this->get_double(name, flag);
  std::string input_units = this->get_string(name + "_units");

  try {
    return units::convert(m_impl->unit_system, value, input_units, units);
  } catch (RuntimeError &e) {
    e.add_context("converting \"%s\" from \"%s\" to \"%s\"",
                  name.c_str(), input_units.c_str(), units.c_str());
    throw;
  }
}

void Config::set_double(const std::string &name, double value,
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

  this->set_double_impl(name, value);
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

Config::Booleans Config::all_booleans() const {
  return this->all_booleans_impl();
}

bool Config::get_boolean(const std::string& name, UseFlag flag) const {
  if (flag == REMEMBER_THIS_USE) {
    m_impl->parameters_used.insert(name);
  }
  return this->get_boolean_impl(name);
}

void Config::set_boolean(const std::string& name, bool value,
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

  this->set_boolean_impl(name, value);
}

static bool special_parameter(const std::string &name) {
  return (ends_with(name, "_doc") or
          ends_with(name, "_units") or
          ends_with(name, "_type") or
          ends_with(name, "_option") or
          ends_with(name, "_choices"));
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

    if (strings[name + "_type"] == "keyword") {
      log.message(v, "  %s%s = \"%s\" (allowed choices: %s)\n",
                  name.c_str(), padding.c_str(), value.c_str(),
                  strings[name + "_choices"].c_str());
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
    double      value = d.second;

    std::string units = strings[name + "_units"]; // will be empty if not set
    std::string padding(max_name_size - name.size(), ' ');

    if (fabs(value) >= 1.0e7 or fabs(value) <= 1.0e-4) {
      // use scientific notation if a number is big or small
      log.message(v, "  %s%s = %13.3e (%s)\n", name.c_str(), padding.c_str(), value, units.c_str());
    } else {
      log.message(v, "  %s%s = %13.5f (%s)\n", name.c_str(), padding.c_str(), value, units.c_str());
    }
  }

  log.message(v,
             "### Booleans:\n"
             "###\n");

  // find max. name size
  max_name_size = 0;
  for (auto b : config.all_booleans()) {
    max_name_size = std::max(max_name_size, b.first.size());
  }

  // print booleans
  for (auto b : config.all_booleans()) {
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

    if (ends_with(p, "_doc")) {
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
 * When called as `set_boolean_from_option(config, "foo", "bar")`,
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
void set_boolean_from_option(Config &config, const std::string &option,
                             const std::string &parameter_name) {

  // get the default value
  bool value = config.get_boolean(parameter_name, Config::FORGET_THIS_USE);
  std::string doc = config.get_string(parameter_name + "_doc", Config::FORGET_THIS_USE);

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

  config.set_boolean(parameter_name, value, CONFIG_USER);
}

//! Sets a configuration parameter from a command-line option.
/*!
  If called as scalar_from_option("foo", "foo"), checks -foo and calls set("foo", value).

  Does nothing if -foo was not set.

  Note that no unit conversion is performed; parameters should be stored in
  input units and converted as needed. (This allows saving parameters without
  converting again.)
*/
void set_scalar_from_option(Config &config, const std::string &name, const std::string &parameter) {
  options::Real option("-" + name,
                       config.get_string(parameter + "_doc", Config::FORGET_THIS_USE),
                       config.get_double(parameter, Config::FORGET_THIS_USE));
  if (option.is_set()) {
    config.set_double(parameter, option, CONFIG_USER);
  }
}

void set_integer_from_option(Config &config, const std::string &name, const std::string &parameter) {
  options::Integer option("-" + name,
                          config.get_string(parameter + "_doc", Config::FORGET_THIS_USE),
                          config.get_double(parameter, Config::FORGET_THIS_USE));
  if (option.is_set()) {
    config.set_double(parameter, option, CONFIG_USER);
  }
}

void set_string_from_option(Config &config, const std::string &name, const std::string &parameter) {

  options::String value("-" + name,
                        config.get_string(parameter + "_doc", Config::FORGET_THIS_USE),
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

  options::Keyword keyword("-" + name,
                           config.get_string(parameter + "_doc", Config::FORGET_THIS_USE),
                           choices,
                           config.get_string(parameter, Config::FORGET_THIS_USE));

  if (keyword.is_set()) {
    config.set_string(parameter, keyword, CONFIG_USER);
  }
}

void set_parameter_from_options(Config &config, const std::string &name) {

  // skip special parameters ("attributes" of parameters)
  if (special_parameter(name)) {
    return;
  }

  // Use parameter name as its own command-line option by default. parameter_name_option can specify
  // a different (possibly shorter) command-line option.
  std::string option = name;

  if (config.is_set(name + "_option")) { // there is a short version of the command-line option
    std::string
      short_option = config.get_string(name + "_option"),
      description  = config.get_string(name + "_doc");

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

  std::string type = "string";
  if (config.is_set(name + "_type")) {
    // will get marked as "used", but that's OK
    type = config.get_string(name + "_type");
  }

  if (type == "string") {
    set_string_from_option(config, option, name);
  } else if (type == "boolean") {
    set_boolean_from_option(config, option, name);
  } else if (type == "scalar") {
    set_scalar_from_option(config, option, name);
  } else if (type == "integer") {
    set_integer_from_option(config, option, name);
  } else if (type == "keyword") {
    // will be marked as "used" and will fail if not set
    std::string choices = config.get_string(name + "_choices");

    set_keyword_from_option(config, option, name, choices);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "parameter type \"%s\" is invalid", type.c_str());
  }
}

void set_config_from_options(Config &config) {
  for (auto d : config.all_doubles()) {
    set_parameter_from_options(config, d.first);
  }

  for (auto s : config.all_strings()) {
    set_parameter_from_options(config, s.first);
  }

  for (auto b : config.all_booleans()) {
    set_parameter_from_options(config, b.first);
  }

  // Energy modeling
  {
    options::Keyword energy("-energy",
                            "choose the energy model (one of 'none', 'cold', 'enthalpy')",
                            "none,cold,enthalpy", "enthalpy");

    if (energy.is_set()) {
      if (energy == "none") {
        config.set_boolean("energy.enabled", false, CONFIG_USER);
        // Allow selecting cold ice flow laws in isothermal mode.
        config.set_boolean("energy.temperature_based", true, CONFIG_USER);
      } else if (energy == "cold") {
        config.set_boolean("energy.enabled", true, CONFIG_USER);
        config.set_boolean("energy.temperature_based", true, CONFIG_USER);
      } else if (energy == "enthalpy") {
        config.set_boolean("energy.enabled", true, CONFIG_USER);
        config.set_boolean("energy.temperature_based", false, CONFIG_USER);
      } else {
        throw RuntimeError(PISM_ERROR_LOCATION, "this can't happen: options::Keyword validates input");
      }
    }
  }

  // -topg_to_phi
  {
    std::vector<double> defaults = {
      config.get_double("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min"),
      config.get_double("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max"),
      config.get_double("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min"),
      config.get_double("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max")
    };

    options::RealList topg_to_phi("-topg_to_phi", "phi_min, phi_max, topg_min, topg_max",
                                  defaults);
    if (topg_to_phi.is_set()) {
      if (topg_to_phi->size() != 4) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "option -topg_to_phi expected 4 numbers; got %d",
                                      (int)topg_to_phi->size());
      }
      config.set_boolean("basal_yield_stress.mohr_coulomb.topg_to_phi.enabled", true);
      config.set_double("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min", topg_to_phi[0]);
      config.set_double("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max", topg_to_phi[1]);
      config.set_double("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min", topg_to_phi[2]);
      config.set_double("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max", topg_to_phi[3]);
    }
  }
  // Ice shelves

  bool nu_bedrock = options::Bool("-nu_bedrock", "constant viscosity near margins");
  if (nu_bedrock) {
    config.set_boolean("stress_balance.ssa.fd.lateral_drag.enabled", true, CONFIG_USER);
  }

  // Shortcuts

  // option "-pik" turns on a suite of PISMPIK effects (but NOT a calving choice,
  // and in particular NOT  "-calving eigen_calving")
  bool pik = options::Bool("-pik", "enable suite of PISM-PIK mechanisms");
  if (pik) {
    config.set_boolean("stress_balance.calving_front_stress_bc", true, CONFIG_USER);
    config.set_boolean("geometry.part_grid.enabled", true, CONFIG_USER);
    config.set_boolean("geometry.remove_icebergs", true, CONFIG_USER);
    config.set_boolean("geometry.grounded_cell_fraction", true, CONFIG_USER);
  }

  if (config.get_string("calving.methods").find("eigen_calving") != std::string::npos) {
    config.set_boolean("geometry.part_grid.enabled", true, CONFIG_USER);
    // eigen-calving requires a wider stencil:
    config.set_double("grid.max_stencil_width", 3);
  }

  // all calving mechanisms require iceberg removal
  if (not config.get_string("calving.methods").empty()) {
    config.set_boolean("geometry.remove_icebergs", true, CONFIG_USER);
  }

  // geometry.remove_icebergs requires part_grid
  if (config.get_boolean("geometry.remove_icebergs")) {
    config.set_boolean("geometry.part_grid.enabled", true, CONFIG_USER);
  }

  bool test_climate_models = options::Bool("-test_climate_models",
                                           "Disable ice dynamics to test climate models");
  if (test_climate_models) {
    config.set_string("stress_balance.model", "none", CONFIG_USER);
    config.set_boolean("energy.enabled", false, CONFIG_USER);
    config.set_boolean("age.enabled", false, CONFIG_USER);
    // let the user decide if they want to use "-no_mass" or not
  }

  // old options
  options::deprecated("-sliding_scale_brutal",
                      "-brutal_sliding' and '-brutal_sliding_scale");
  options::deprecated("-ssa_sliding", "-stress_balance ...");
  options::deprecated("-ssa_floating_only", "-stress_balance ...");
  options::deprecated("-sia", "-stress_balance ...");
  options::deprecated("-no_sia", "-stress_balance ...");
  options::deprecated("-hold_tauc", "-yield_stress constant");
  options::deprecated("-ocean_kill", "-calving ocean_kill -ocean_kill_file foo.nc");
  options::deprecated("-eigen_calving", "-calving eigen_calving -eigen_calving_K XXX");
  options::deprecated("-calving_at_thickness",
                      "-calving thickness_calving -thickness_calving_threshold XXX");
  options::deprecated("-float_kill", "-calving float_kill");
  options::deprecated("-no_energy", "-energy none");
  options::deprecated("-cold", "-energy cold");
  options::deprecated("-boot_file", "-bootstrap -i");
}

//! Create a configuration database using command-line options.
Config::Ptr config_from_options(MPI_Comm com, const Logger &log, units::System::Ptr sys) {

  DefaultConfig::Ptr config(new DefaultConfig(com, "pism_config", "-config", sys)),
    overrides(new DefaultConfig(com, "pism_overrides", "-config_override", sys));
  overrides->init(log);
  config->init_with_default(log);
  config->import_from(*overrides);
  set_config_from_options(*config);

  return config;
}

ConfigWithPrefix::ConfigWithPrefix(Config::ConstPtr c, const std::string &prefix)
  : m_prefix(prefix), m_config(c) {
  // empty
}

double ConfigWithPrefix::get_double(const std::string &name) const {
  return m_config->get_double(m_prefix + name);
}

double ConfigWithPrefix::get_double(const std::string &name, const std::string &units) const {
  return m_config->get_double(m_prefix + name, units);
}

std::string ConfigWithPrefix::get_string(const std::string &name) const {
  return m_config->get_string(m_prefix + name);
}

bool ConfigWithPrefix::get_boolean(const std::string& name) const {
  return m_config->get_boolean(m_prefix + name);
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

  for (auto p : all_booleans()) {
    result.insert(p.first);
  }

  return result;
}

} // end of namespace pism
