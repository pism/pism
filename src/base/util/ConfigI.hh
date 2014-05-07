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

#ifndef _CONFIGI_H_
#define _CONFIGI_H_

#include <string>

namespace pism {

class ConfigI {
public:
  virtual ~ConfigI() {}

  virtual void set_double(std::string const &name, double value) = 0;
  virtual void set_string(std::string const &name, std::string const &value) = 0;
  virtual void   set_flag(std::string const &, bool) = 0;

  virtual double get(std::string const &) const = 0;
  inline double get_double(std::string const &name) const
    { return get(name); }
  virtual bool   get_flag(std::string const &) const = 0;
  virtual std::string get_string(std::string const &name) const = 0;

};

}  // namespace pism

#endif /* _CONFIGI_H_ */
