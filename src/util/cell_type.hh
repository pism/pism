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
  UNKNOWN        = -1,
  ICE_FREE_LAND = 0,
  ICY_LAND      = 1,
  // values associated with "land" are less than ones for "water"
  ICE_FREE_LAKE  = 2,
  ICY_LAKE       = 3,
  ICE_FREE_OCEAN = 4,
  ICY_OCEAN      = 5
};

//! An wet cell (floating ice or ice-free).
inline bool water(int M) {
  return M > ICY_LAND;
}

//! Grounded cell (grounded ice or ice-free).
inline bool land(int M) {
  return M <= ICY_LAND;
}

//! Ice-filled cell (grounded or floating).
inline bool icy(int M) {
  return (M == ICY_LAND) or (M == ICY_LAKE) or (M == ICY_OCEAN);
}

inline bool grounded_ice(int M) {
  return M == ICY_LAND;
}

inline bool floating_ice(int M) {
  return (M == ICY_LAKE) or (M == ICY_OCEAN);
}

//! Ice-free cell (grounded or ocean).
inline bool ice_free(int M) {
  return (M == ICE_FREE_LAND) or (M == ICE_FREE_LAKE) or (M == ICE_FREE_OCEAN);
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

inline bool lake(int M) {
  return (M == ICE_FREE_LAKE) or (M == ICY_LAKE);
}

inline bool ocean(int M) {
  return (M == ICE_FREE_OCEAN) or (M == ICY_OCEAN);
}
} // namespace cell_type

} // end of namespace pism

#endif /* PISM_CELL_TYPE_H */
