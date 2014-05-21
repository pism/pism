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

#ifndef _CONFIGJSON_H_
#define _CONFIGJSON_H_

#include "Config.hh"
#include <jansson.h>
#include <vector>

namespace pism {

/*! @brief The JSON-based Config implementation.
 *
 */
class ConfigJSON : public Config {
public:
  ConfigJSON();
  virtual ~ConfigJSON();
  int init_from_file(const std::string &filename);
  int init_from_string(const std::string &string);
  std::string dump() const;
private:
  virtual void set_double_impl(const std::string &name, double value);
  virtual void set_boolean_impl(const std::string &name, bool value);
  virtual void set_string_impl(const std::string &name, const std::string &value);

  virtual double get_double_impl(const std::string &name) const;
  virtual std::string get_string_impl(const std::string &name) const;
  virtual bool get_boolean_impl(const std::string &name) const;
  json_t *m_data;
};

} // end of namespace pism

#endif /* _CONFIGJSON_H_ */
