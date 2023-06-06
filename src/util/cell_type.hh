// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2022, 2023 Constantine Khroulev and David Maxwell
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

const int CELL_TYPE_ICY   = 1;
const int CELL_TYPE_LAND  = (1 << 1);
const int CELL_TYPE_LAKE  = (1 << 2);
const int CELL_TYPE_OCEAN = (1 << 3);

enum Value : int {
  UNKNOWN        = 0,
  ICE_FREE_LAND  = CELL_TYPE_LAND,
  ICY_LAND       = CELL_TYPE_LAND | CELL_TYPE_ICY,
  ICE_FREE_LAKE  = CELL_TYPE_LAKE,
  ICY_LAKE       = CELL_TYPE_LAKE | CELL_TYPE_ICY,
  ICE_FREE_OCEAN = CELL_TYPE_OCEAN,
  ICY_OCEAN      = CELL_TYPE_OCEAN | CELL_TYPE_ICY
};

//! \brief An wet cell (floating ice or ice-free).
  inline bool wet(int M) {
    return ((M & CELL_TYPE_LAKE) != 0) or ((M & CELL_TYPE_OCEAN) != 0);
  }
  //! \brief Grounded cell (grounded ice or ice-free).
  inline bool grounded(int M) {
    return (M & CELL_TYPE_LAND) != 0;
  }
  //! \brief Ice-filled cell (grounded or floating).
  inline bool icy(int M) {
    return (M & CELL_TYPE_ICY) != 0;
  }
  inline bool grounded_ice(int M) {
    return M == ICY_LAND;
  }
  inline bool floating_ice(int M) {
    return (M == ICY_LAKE) or (M == ICY_OCEAN);
  }
  //! \brief Ice-free cell (grounded or ocean).
  inline bool ice_free(int M) {
    return not icy(M);
  }
  inline bool ice_free_ocean(int M) {
    return M == ICE_FREE_OCEAN;
  }
  inline bool ice_free_land(int M) {
    return M == ICE_FREE_LAND;
  }
}

} // end of namespace pism

#endif /* PISM_CELL_TYPE_H */
