/* Copyright (C) 2016 PISM Authors
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

namespace pism {

CalvingFrontRetreat::CalvingFrontRetreat(IceGrid::ConstPtr g)
  : Component(g) {

}

CalvingFrontRetreat::~CalvingFrontRetreat() {

}

//! Update ice thickness, ice volume in partially-filled cells (Href),
//! and cell type mask by applying the 2D horizontal calving rate.
void CalvingFrontRetreat::update(double dt,
                                 IceModelVec2CellType &pism_mask,
                                 IceModelVec2S &Href,
                                 IceModelVec2S &ice_thickness) {

}


} // end of namespace pism
