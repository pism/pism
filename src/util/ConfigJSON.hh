/* Copyright (C) 2014, 2016, 2017, 2018, 2025 PISM Authors
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

#ifndef _CONFIGJSON_H_
#define _CONFIGJSON_H_

#include <jansson.h>

#include "pism/util/Config.hh"

namespace pism {

/*! @brief The JSON-based Config implementation.
 *
 */
class ConfigJSON : public Config {
public:
  ConfigJSON(units::System::Ptr unit_system);
  virtual ~ConfigJSON();

  void init_from_file(const std::string &filename);
  void init_from_string(const std::string &string);
  std::string dump() const;

private:
  void read_impl(const File &nc);
  void write_impl(const File &nc) const;

  bool is_set_impl(const std::string &name) const;

  Doubles all_doubles_impl() const;
  double get_number_impl(const std::string &name) const;
  std::vector<double> get_numbers_impl(const std::string &name) const;

  void set_number_impl(const std::string &name, double value);
  void set_numbers_impl(const std::string &name,
                                const std::vector<double> &values);

  Strings all_strings_impl() const;
  std::string get_string_impl(const std::string &name) const;
  void set_string_impl(const std::string &name, const std::string &value);

  Flags all_flags_impl() const;

  bool get_flag_impl(const std::string& name) const;
  void set_flag_impl(const std::string& name, bool value);

  json_t *m_data;
};

} // end of namespace pism

#endif /* _CONFIGJSON_H_ */
