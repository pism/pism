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

#ifndef _OPTIONS_H_
#define _OPTIONS_H_

namespace pism {
namespace options {

template <typename T>
class Option {
public:
  Option() {
    m_is_set = false;
  }
  operator T() {
    return m_value;
  }
  bool is_set() {
    return m_is_set;
  }
  T value() {
    return m_value;
  }
protected:
  T m_value;
  bool m_is_set;
  void set(T value, bool is_set) {
    m_value = value;
    m_is_set = is_set;
  }
};

} // end of namespace options
} // end of namespace pism


#endif /* _OPTIONS_H_ */
