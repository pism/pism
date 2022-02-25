/* Copyright (C) 2016, 2017, 2018, 2019, 2021, 2022 PISM Authors
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

#ifndef FRONTRETREAT_H
#define FRONTRETREAT_H

#include "pism/util/Component.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {

class Geometry;

//! An abstract class implementing calving front retreat resulting from application of a
//! spatially-variable horizontal retreat rate.
/*! The retreat rate may correspond to frontal melting or calving. Requires the
    "part_grid" mechanism and a grid with dx == dy.
 */
class FrontRetreat : public Component {
public:
  FrontRetreat(IceGrid::ConstPtr g);
  ~FrontRetreat() = default;

  void update_geometry(double dt,
                       const Geometry &geometry,
                       const IceModelVec2S &bc_mask,
                       const IceModelVec2S &retreat_rate,
                       IceModelVec2S &Href,
                       IceModelVec2S &ice_thickness);

  MaxTimestep max_timestep(const CellTypeArray1 &cell_type,
                           const IceModelVec2S &bc_mask,
                           const IceModelVec2S &retreat_rate) const;
private:

  void compute_modified_mask(const CellTypeArray1 &input,
                             CellTypeArray1 &output) const;

  // Ghosted cell type mask
  //
  // We make a copy here because frontal retreat code uses a modified mask if
  // geometry.front_retreat.wrap_around is false.
  CellTypeArray1 m_cell_type;
  // Temporary storage for distributing ice loss to "full" (as opposed to "partially
  // filled") cells near the front
  Array2SGhosted<1> m_tmp;
};

} // end of namespace pism


#endif /* FRONTRETREAT_H */
