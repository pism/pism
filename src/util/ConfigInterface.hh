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

#ifndef _PISMCONFIGINTERFACE_H_
#define _PISMCONFIGINTERFACE_H_

#include <memory>
#include <set>
#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include "Units.hh"

namespace pism {

class PIO;
class Logger;


//! Flag used by `set_...()` methods.
/** Meanings:
 *
 * - `DEFAULT`: set the default value; has no effect if a parameter was set by a user at the time
 *   of the call
 * - `FORCE`: forcibly set a parameter; unconditionally overrides previous values
 * - `USER`: forcibly set a parameter; unconditionally overrides previous values and marks this
 *   parameter as set by the user. This affects future `set_...()` calls using the `DEFAULT` flag
 *   value and results of `parameters_set_by_user()`.
 */
enum ConfigSettingFlag {CONFIG_DEFAULT = 0, CONFIG_FORCE = 1, CONFIG_USER = 2};


//! A class for storing and accessing PISM configuration flags and parameters.
class Config {
public:
  typedef std::shared_ptr<Config> Ptr;
  typedef std::shared_ptr<const Config> ConstPtr;

  Config(units::System::Ptr unit_system);
  virtual ~Config();

  //! Flag used by `get_...()` methods.
  /** Meanings:
   *
   * - `REMEMBER_THIS_USE` (the default): add the name of a parameter to the list of parameters used
   *    by a model run.
   * - `FORGET_THIS_USE`: don't add the name of a parameter to the list of used parameters. This is
   *    necessary to be able to get the current value of a parameter to be used as the default when
   *    processing a command-line option.
   */
  enum UseFlag {REMEMBER_THIS_USE = 0, FORGET_THIS_USE = 1};

  // Import settings from an override file
  void import_from(const Config &other);

  const std::set<std::string>& parameters_set_by_user() const;
  const std::set<std::string>& parameters_used() const;

  void read(MPI_Comm com, const std::string &filename);
  void write(MPI_Comm com, const std::string &filename, bool append = true) const;
  std::string filename() const;

  void read(const PIO &nc);
  void write(const PIO &nc) const;

  bool is_set(const std::string &name) const;

  // doubles
  typedef std::map<std::string, double> Doubles;
  Doubles all_doubles() const;

  double get_double(const std::string &name, UseFlag flag = REMEMBER_THIS_USE) const;
  double get_double(const std::string &name, const std::string &units,
                    UseFlag flag = REMEMBER_THIS_USE) const;
  void set_double(const std::string &name, double value, ConfigSettingFlag flag = CONFIG_FORCE);

  // strings
  typedef std::map<std::string, std::string> Strings;
  Strings all_strings() const;

  std::string get_string(const std::string &name, UseFlag flag = REMEMBER_THIS_USE) const;
  void set_string(const std::string &name, const std::string &value, ConfigSettingFlag flag = CONFIG_FORCE);

  // booleans
  typedef std::map<std::string, bool> Booleans;
  Booleans all_booleans() const;

  std::set<std::string> keys() const;

  bool get_boolean(const std::string& name, UseFlag flag = REMEMBER_THIS_USE) const;
  void set_boolean(const std::string& name, bool value, ConfigSettingFlag flag = CONFIG_FORCE);

  // Implementations
protected:
  virtual void read_impl(const PIO &nc) = 0;
  virtual void write_impl(const PIO &nc) const = 0;

  virtual bool is_set_impl(const std::string &name) const = 0;

  virtual Doubles all_doubles_impl() const = 0;
  virtual double get_double_impl(const std::string &name) const = 0;
  virtual void set_double_impl(const std::string &name, double value) = 0;

  virtual Strings all_strings_impl() const = 0;
  virtual std::string get_string_impl(const std::string &name) const = 0;
  virtual void set_string_impl(const std::string &name, const std::string &value) = 0;

  virtual Booleans all_booleans_impl() const = 0;

  virtual bool get_boolean_impl(const std::string& name) const = 0;
  virtual void set_boolean_impl(const std::string& name, bool value) = 0;
private:
  struct Impl;
  Impl *m_impl;
};

class ConfigWithPrefix {
public:
  ConfigWithPrefix(Config::ConstPtr c, const std::string &prefix);

  double get_double(const std::string &name) const;
  double get_double(const std::string &name, const std::string &units) const;

  std::string get_string(const std::string &name) const;

  bool get_boolean(const std::string& name) const;

  void reset_prefix(const std::string &prefix);

private:
  std::string m_prefix;
  Config::ConstPtr m_config;
};

Config::Ptr config_from_options(MPI_Comm com, const Logger &log, units::System::Ptr unit_system);

//! Set configuration parameters using command-line options.
void set_config_from_options(Config &config);

//! Set one parameter using command-line options.
void set_parameter_from_options(Config &config, const std::string &name);

//! Set one boolean parameter using command-line options.
void set_boolean_from_option(Config &config,
                             const std::string &option,const std::string &parameter);

//! Set one scalar parameter using command-line options.
void set_scalar_from_option(Config &config,
                            const std::string &option, const std::string &parameter);

//! Set one free-form string parameter using command-line options.
void set_string_from_option(Config &config,
                            const std::string &option, const std::string &parameter);

//! Set one keyword parameter using command-line options.
void set_keyword_from_option(Config &config,
                             const std::string &option, const std::string &parameter,
                             const std::string &choices);

//! Report configuration parameters to `stdout`.
void print_config(const Logger &log, int verbosity_threshhold, const Config &config);

//! Report unused configuration parameters to `stdout`.
void print_unused_parameters(const Logger &log, int verbosity_threshhold,
                             const Config &config);

} // end of namespace pism

#endif /* _PISMCONFIGINTERFACE_H_ */
