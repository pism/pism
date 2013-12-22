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
#include <string>
#include <set>

//! A class for reading, writing and accessing PISM configuration flags and parameters.
class PISMConfig {
public:
  PISMConfig(MPI_Comm com, std::string name, PISMUnitSystem unit_system);
  ~PISMConfig();

  void set_double(std::string name, double value);
  void set_string(std::string name, std::string value);

  bool is_set(std::string name) const;

  PetscErrorCode print_to_stdout(PetscInt verbosity_threshhold = 4) const;
  PetscErrorCode warn_about_unused_parameters() const;
  PetscErrorCode read(const PIO &nc);
  PetscErrorCode write(const PIO &nc) const;

  PetscErrorCode read(std::string filename);
  PetscErrorCode write(std::string filename, bool append = true) const;

  std::string get_config_filename() const;
  PISMUnitSystem get_unit_system() const;
  double get(std::string) const;
  double get(std::string name, std::string u1, std::string u2) const;
  bool   get_flag(std::string) const;
  std::string get_string(std::string name) const;
  // Set a flag (overriding the default in pism_config.cdl). Should not be used
  // in pismr code.
  void   set_flag(std::string, bool);
  // Set parameters and remember that they were set using a command-line option
  PetscErrorCode set_flag_from_option(std::string name, bool value);
  PetscErrorCode set_scalar_from_option(std::string name, double value);
  PetscErrorCode set_string_from_option(std::string name, std::string value);
  PetscErrorCode set_keyword_from_option(std::string name, std::string value);
  // Set parameters by ptocessing a command-line option
  PetscErrorCode flag_from_option(std::string, std::string);
  PetscErrorCode scalar_from_option(std::string, std::string);
  PetscErrorCode string_from_option(std::string, std::string);
  PetscErrorCode keyword_from_option(std::string, std::string, std::string);
  // Import settings from an override file
  void import_from(const PISMConfig &other);
  void update_from(const PISMConfig &other);

private:
  MPI_Comm m_com;
  PISMUnitSystem m_unit_system;
  NCVariable m_data;
  std::string m_config_filename;
  //!< \brief the name of the file this config database was initialized from 
  double get_quiet(std::string name) const;
  std::string get_string_quiet(std::string name) const;
  bool   get_flag_quiet(std::string name) const;
  const NCVariable& get_data() const;

  std::set<std::string> m_parameters_set;
  mutable std::set<std::string> m_parameters_used;
  bool m_options_left_set;
};

#endif /* _PISMCONFIG_H_ */
