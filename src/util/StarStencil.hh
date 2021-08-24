/* Copyright (C) 2015, 2016, 2021 PISM Authors
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

#ifndef _STARSTENCIL_H_
#define _STARSTENCIL_H_

namespace pism {
enum Direction {North = 0, East, South, West};

//! \brief Star stencil points (in the map-plane).
template <typename T>
struct StarStencil {
  T ij, e, w, n, s;

  StarStencil() = default;

  StarStencil(T value) {
    set(value);
  }

  void set(T input) {
    ij = e = w = n = s = input;
  }

  //! Get the element corresponding to a given direction.
  //! Use foo.ij to get the value at i,j (center of the star).
  inline T& operator[](Direction direction) {
    switch (direction) {
    default:                    // just to silence the warning
    case North:
      return n;
    case East:
      return e;
    case South:
      return s;
    case West:
      return w;
    }
  }

  inline const T& operator[](Direction direction) const {
    switch (direction) {
    default:                    // just to silence the warning
    case North:
      return n;
    case East:
      return e;
    case South:
      return s;
    case West:
      return w;
    }
  }
};

}

#endif /* _STARSTENCIL_H_ */
