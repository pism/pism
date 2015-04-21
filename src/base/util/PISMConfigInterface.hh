/* Copyright (C) 2015 PISM Authors
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

#include <set>
#include <map>
#include <string>
#include <mpi.h>

#include "PISMUnits.hh"
#include "pism_memory.hh"

namespace pism {

class PIO;

//! A class for storing and accessing PISM configuration flags and parameters.
class Config {
public:
  typedef PISM_SHARED_PTR_NSPACE::shared_ptr<Config> Ptr;
  typedef PISM_SHARED_PTR_NSPACE::shared_ptr<const Config> ConstPtr;

  Config(units::System::Ptr unit_system);
  virtual ~Config();

  enum SetBy {SYSTEM = 0, USER = 1};

  // methods implemented in the base class
  units::System::Ptr unit_system() const;

  // Import settings from an override file
  void import_from(const Config &other);
  // Update settings stored in an "override" config to make sure that when we save an override
  // config it contains values that were used by a run.
  void update_from(const Config &other);

  const std::set<std::string>& parameters_set_by_user() const;
  const std::set<std::string>& parameters_used() const;

  void read(MPI_Comm com, const std::string &filename);
  void write(MPI_Comm com, const std::string &filename, bool append = true) const;
  std::string filename() const;

  // end of methods implemented in the base class

  void read(const PIO &nc);
  void write(const PIO &nc) const;

  bool is_set(const std::string &name) const;

  // doubles
  typedef std::map<std::string, double> Doubles;
  Doubles all_doubles() const;

  double get_double(const std::string &name) const;
  double get_double(const std::string &name, const std::string &u1, const std::string &u2) const;
  void set_double(const std::string &name, double value, SetBy flag = SYSTEM);

  // strings
  typedef std::map<std::string, std::string> Strings;
  Strings all_strings() const;

  std::string get_string(const std::string &name) const;
  void set_string(const std::string &name, const std::string &value, SetBy flag = SYSTEM);

  // booleans
  typedef std::map<std::string, bool> Booleans;
  Booleans all_booleans() const;

  bool get_boolean(const std::string& name) const;
  void set_boolean(const std::string& name, bool value, SetBy flag = SYSTEM);

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
protected:
  // Set parameters by ptocessing a command-line option
  void boolean_from_option(const std::string &option, const std::string &parameter);
  void scalar_from_option(const std::string &option, const std::string &parameter);
  void string_from_option(const std::string &option, const std::string &parameter);
  void keyword_from_option(const std::string &option, const std::string &parameter,
                           const std::string &choices);

  void set_from_options();
private:
  struct Impl;
  Impl *m_impl;
};

void print_config(int verbosity_threshhold, MPI_Comm com, const Config &config);
void print_unused_parameters(int verbosity_threshhold, MPI_Comm com, const Config &config);

} // end of namespace pism

#endif /* _PISMCONFIGINTERFACE_H_ */
