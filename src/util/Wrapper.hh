/* Copyright (C) 2015, 2016, 2021, 2025 PISM Authors
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

#ifndef _WRAPPER_H_
#define _WRAPPER_H_

namespace pism {

template<typename T>
class Wrapper {
public:
  operator T() const {
    return m_value;
  }
  T get() const {
    return m_value;
  }
  T* rawptr() {
    return &m_value;
  }
protected:
  Wrapper() {
    // empty
  }
  T m_value;
private:
  Wrapper(Wrapper const &);
  Wrapper & operator=(Wrapper const &);
};

} // end of namespace pism

#endif /* _WRAPPER_H_ */
