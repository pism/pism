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

#ifndef _PISMCONFIG_H_
#define _PISMCONFIG_H_

#include "NCVariable.hh"
#include "PISMUnits.hh"
#include "ConfigI.hh"
#include <string>
#include <set>

namespace pism {

//! A class for reading, writing and accessing PISM configuration flags and parameters.
class Config : public ConfigI {
public:
  Config(MPI_Comm com, const std::string &name, const UnitSystem &unit_system);
  ~Config();

  void set_double(const std::string &name, double value);
  void set_string(const std::string &name, const std::string &value);

  bool is_set(const std::string &name) const;

  PetscErrorCode print_to_stdout(int verbosity_threshhold = 4) const;
  PetscErrorCode warn_about_unused_parameters() const;
  PetscErrorCode read(const PIO &nc);
  PetscErrorCode write(const PIO &nc) const;

  PetscErrorCode read(const std::string &filename);
  PetscErrorCode write(const std::string &filename, bool append = true) const;

  std::string get_config_filename() const;
  UnitSystem get_unit_system() const;
  double get(const std::string &) const;
  double get(const std::string &name, const std::string &u1, const std::string &u2) const;
  bool   get_flag(const std::string&) const;
  std::string get_string(const std::string &name) const;
  // Set a flag (overriding the default in pism_config.cdl). Should not be used
  // in pismr code.
  void   set_flag(const std::string&, bool);
  // Set parameters and remember that they were set using a command-line option
  PetscErrorCode set_flag_from_option(const std::string &name, bool value);
  PetscErrorCode set_scalar_from_option(const std::string &name, double value);
  PetscErrorCode set_string_from_option(const std::string &name, const std::string &value);
  PetscErrorCode set_keyword_from_option(const std::string &name, const std::string &value);
  // Set parameters by ptocessing a command-line option
  PetscErrorCode flag_from_option(const std::string &, const std::string &);
  PetscErrorCode scalar_from_option(const std::string &, const std::string &);
  PetscErrorCode string_from_option(const std::string &, const std::string &);
  PetscErrorCode keyword_from_option(const std::string &, const std::string &, const std::string &);
  // Import settings from an override file
  void import_from(const Config &other);
  void update_from(const Config &other);

private:
  MPI_Comm m_com;
  UnitSystem m_unit_system;
  NCVariable m_data;
  std::string m_config_filename;
  //!< \brief the name of the file this config database was initialized from 
  double get_quiet(const std::string &name) const;
  std::string get_string_quiet(const std::string &name) const;
  bool   get_flag_quiet(const std::string &name) const;
  const NCVariable& get_data() const;

  std::set<std::string> m_parameters_set;
  mutable std::set<std::string> m_parameters_used;
  bool m_options_left_set;
};

} // end of namespace pism

#endif /* _PISMCONFIG_H_ */
