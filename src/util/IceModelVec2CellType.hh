/* Copyright (C) 2016, 2019 PISM Authors
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

#ifndef ICEMODELVEC2CELLTYPE_H
#define ICEMODELVEC2CELLTYPE_H

#include "iceModelVec.hh"
#include "Mask.hh"

namespace pism {

//! "Cell type" mask. Adds convenience methods to IceModelVec2Int.
class IceModelVec2CellType : public IceModelVec2Int {
public:

  typedef std::shared_ptr<IceModelVec2CellType> Ptr;
  typedef std::shared_ptr<const IceModelVec2CellType> ConstPtr;
  IceModelVec2CellType()
    : IceModelVec2Int() {
    // empty
  }

  IceModelVec2CellType(IceGrid::ConstPtr grid, const std::string &name,
                       IceModelVecKind ghostedp, int width = 1)
    : IceModelVec2Int(grid, name, ghostedp, width) {
    // empty
  }

  inline bool ocean(int i, int j) const {
    return mask::ocean(as_int(i, j));
  }

  inline bool grounded(int i, int j) const {
    return mask::grounded(as_int(i, j));
  }

  inline bool icy(int i, int j) const {
    return mask::icy(as_int(i, j));
  }

  inline bool grounded_ice(int i, int j) const {
    return mask::grounded_ice(as_int(i, j));
  }

  inline bool floating_ice(int i, int j) const {
    return mask::floating_ice(as_int(i, j));
  }

  inline bool ice_free(int i, int j) const {
    return mask::ice_free(as_int(i, j));
  }

  inline bool ice_free_ocean(int i, int j) const {
    return mask::ice_free_ocean(as_int(i, j));
  }

  inline bool ice_free_land(int i, int j) const {
    return mask::ice_free_land(as_int(i, j));
  }

  inline bool ice_free_enclosed_ocean(int i, int j) const {
    return mask::ice_free_enclosed_ocean(as_int(i, j));
  }

  inline bool ice_free_open_ocean(int i, int j) const {
    return mask::ice_free_open_ocean(as_int(i, j));
  }

  //! \brief Ice margin (ice-filled with at least one of four neighbors ice-free).
  inline bool ice_margin(int i, int j) const {
    return icy(i, j) and (ice_free(i + 1, j) or ice_free(i - 1, j) or
                          ice_free(i, j + 1) or ice_free(i, j - 1));
  }

  //! \brief Ice-free margin (at least one of four neighbors has ice).
  inline bool next_to_ice(int i, int j) const {
    return (icy(i + 1, j) or icy(i - 1, j) or icy(i, j + 1) or icy(i, j - 1));
  }

  inline bool next_to_floating_ice(int i, int j) const {
    return (floating_ice(i + 1, j) or floating_ice(i - 1, j) or
            floating_ice(i, j + 1) or floating_ice(i, j - 1));
  }

  inline bool next_to_grounded_ice(int i, int j) const {
    return (grounded_ice(i + 1, j) or grounded_ice(i - 1, j) or
            grounded_ice(i, j + 1) or grounded_ice(i, j - 1));
  }

  inline bool next_to_ice_free_land(int i, int j) const {
    return (ice_free_land(i + 1, j) or ice_free_land(i - 1, j) or
            ice_free_land(i, j + 1) or ice_free_land(i, j - 1));
  }

  inline bool next_to_ice_free_ocean(int i, int j) const {
    return (ice_free_ocean(i + 1, j) or ice_free_ocean(i - 1, j) or
            ice_free_ocean(i, j + 1) or ice_free_ocean(i, j - 1));
  }

  inline bool next_to_ice_free_open_ocean(int i, int j) const {
    return (ice_free_open_ocean(i + 1, j) or ice_free_open_ocean(i - 1, j) or
            ice_free_open_ocean(i, j + 1) or ice_free_open_ocean(i, j - 1));
  }
};

} // end of namespace pism


#endif /* ICEMODELVEC2CELLTYPE_H */
