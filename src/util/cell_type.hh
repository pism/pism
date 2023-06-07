// Copyright (C) 2011--2023 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef PISM_CELL_TYPE_H
#define PISM_CELL_TYPE_H

namespace pism {

namespace array {
class Scalar;
}

namespace cell_type {

enum Value : int {
  UNKNOWN        = -2,          // FIXME: this will be interpreted as "ice free land" by the code below
  ICE_FREE_LAND = 0, // even means "ice free"
  ICY_LAND      = 1, // odd means "icy"
  // values associated with "land" are less than ones for "water"
  ICE_FREE_LAKE  = 2, // even means "ice free"
  ICY_LAKE       = 3, // odd means "icy"
  ICE_FREE_OCEAN = 4, // even means "ice free"
  ICY_OCEAN      = 5  // odd means "icy"
};

//! \brief An wet cell (floating ice or ice-free).
  inline bool water(int M) {
    return M > ICY_LAND;
  }
  //! \brief Grounded cell (grounded ice or ice-free).
  inline bool land(int M) {
    return M <= ICY_LAND;
  }
  //! \brief Ice-filled cell (grounded or floating).
  inline bool icy(int M) {
    return M % 2 == 1;          // odd means "icy"
  }
  inline bool grounded_ice(int M) {
    return M == ICY_LAND;
  }
  inline bool floating_ice(int M) {
    return (M == ICY_LAKE) or (M == ICY_OCEAN);
  }
  //! \brief Ice-free cell (grounded or ocean).
  inline bool ice_free(int M) {
    return M % 2 == 0;
  }
  // inline bool ice_free_ocean(int M) {
  //   return M == ICE_FREE_OCEAN;
  // }
  inline bool ice_free_water(int M) {
    return (M == ICE_FREE_LAKE) or (M == ICE_FREE_OCEAN);
  }
  inline bool ice_free_land(int M) {
    return M == ICE_FREE_LAND;
  }
  } // namespace cell_type

} // end of namespace pism

#endif /* PISM_CELL_TYPE_H */
