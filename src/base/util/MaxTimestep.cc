/* Copyright (C) 2015, 2016 PISM Authors
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

#include "MaxTimestep.hh"

namespace pism {

// Time step restrictions
MaxTimestep::MaxTimestep()
  : m_is_finite(false), m_value(0.0) {
  // empty
}

MaxTimestep::MaxTimestep(double v)
  : m_is_finite(true), m_value(v) {
  // empty
}

bool MaxTimestep::is_finite() const {
  return m_is_finite;
}

double MaxTimestep::value() const {
  return m_value;
}

bool operator==(const MaxTimestep &a, const MaxTimestep &b) {
  return (a.is_finite() == b.is_finite()) and (a.value() == b.value());
}

bool operator<(const MaxTimestep &a, const MaxTimestep &b) {
  if (a.is_finite() and b.is_finite()) {
    return a.value() < b.value();
  } else if (a.is_finite()) {
    return true;
  } else if (b.is_finite()) {
    return false;
  } else {
    return false;
  }
}

bool operator>(const MaxTimestep &a, const MaxTimestep &b) {
  return (not (a == b)) and (not (a < b));
}

} // end of namespace pism
