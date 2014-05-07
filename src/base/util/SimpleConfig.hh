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

#ifndef _SIMPLECONFIG_H_
#define _SIMPLECONFIG_H_

#include <string>
#include <map>
#include <iostream>
#include "ConfigI.hh"

namespace pism {

class SimpleConfig : public ConfigI
{

  friend std::ostream &operator<<(std::ostream &out, SimpleConfig const &config);

public:
  void set_double(std::string const &name, double value)
	{ doubles.insert(std::make_pair(name, value)); }
  void set_string(std::string const &name, std::string const &value)
	{ strings.insert(std::make_pair(name, value)); }
  void   set_flag(std::string const &name, bool value)
	{ bools.insert(std::make_pair(name, value)); }

  double get(std::string const &name) const
	{ return doubles.find(name)->second; }
  std::string get_string(std::string const &name) const
	{ return strings.find(name)->second; }
  bool   get_flag(std::string const &name) const
	{ return bools.find(name)->second; }

private:
  std::map<std::string, double> doubles;
  std::map<std::string, std::string> strings;
  std::map<std::string, bool> bools;
};

extern std::ostream &operator<<(std::ostream &out, SimpleConfig const &config);

}  // namespace pism

#endif /* _SIMPLECONFIG_H_ */
