/* Copyright (C) 2014, 2015, 2016, 2017, 2018, 2021, 2024, 2025 PISM Authors
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

#ifndef _PISM_NETCDF_CONFIG_H_
#define _PISM_NETCDF_CONFIG_H_

#include <string>
#include <set>

#include "pism/util/Config.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

class Logger;

//! A class for reading, writing and accessing PISM configuration flags and parameters.
class NetCDFConfig : public Config {
public:
  NetCDFConfig(const std::string &name, units::System::Ptr unit_system);
  ~NetCDFConfig();

protected:
  void read_impl(const File &file);
  void write_impl(const File &file) const;

  bool is_set_impl(const std::string &name) const;

  // doubles
  Doubles all_doubles_impl() const;
  double get_number_impl(const std::string &name) const;
  std::vector<double> get_numbers_impl(const std::string &name) const;

  void set_number_impl(const std::string &name, double value);
  void set_numbers_impl(const std::string &name,
                        const std::vector<double> &values);
  // strings
  Strings all_strings_impl() const;
  std::string get_string_impl(const std::string &name) const;
  void set_string_impl(const std::string &name, const std::string &value);

  // flags
  Flags all_flags_impl() const;
  bool get_flag_impl(const std::string& name) const;
  void set_flag_impl(const std::string& name, bool value);
private:
  VariableMetadata m_data;
  //! @brief the name of the file this config database was initialized from
  std::string m_config_filename;
};

} // end of namespace pism

#endif /* _PISM_NETCDF_CONFIG_H_ */
