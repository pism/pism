/* Copyright (C) 2018 PISM Authors
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

#include "PicoGeometry.hh"
#include "pism/calving/connected_components.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace ocean {

PicoGeometry::PicoGeometry(IceGrid::ConstPtr grid)
  : Component(grid),
    m_continental_shelf(grid, "pico_ocean_contshelf_mask", WITHOUT_GHOSTS),
    m_boxes(grid, "pico_ocean_box_mask", WITHOUT_GHOSTS),
    m_ice_shelves(grid, "pico_shelf_mask", WITHOUT_GHOSTS),
    m_tmp(grid, "temporary_storage", WITHOUT_GHOSTS) {

  m_tmp_p0 = m_tmp.allocate_proc0_copy();
}

PicoGeometry::~PicoGeometry() {
  // empty
}

const IceModelVec2Int& PicoGeometry::continental_shelf_mask() const {
  return m_continental_shelf;
}

const IceModelVec2Int& PicoGeometry::box_mask() const {
  return m_boxes;
}

const IceModelVec2Int& PicoGeometry::ice_shelf_mask() const {
  return m_ice_shelves;
}

void PicoGeometry::update(const IceModelVec2S &bed_elevation,
                          const IceModelVec2CellType &cell_type) {

}

void PicoGeometry::compute_ice_rises(const IceModelVec2CellType &cell_type,
                                     IceModelVec2Int &result) {
  IceModelVec::AccessList list{&cell_type, &m_tmp};

  // mask of zeros and ones: one if grounded ice, zero otherwise
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.grounded(i, j)) {
      m_tmp(i, j) = 1.0;
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  m_tmp.put_on_proc0(*m_tmp_p0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      petsc::VecArray mask_p0(*m_tmp_p0);
      label_connected_components(mask_p0.get(), m_grid->My(), m_grid->Mx(), false, 0.0);
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_tmp.get_from_proc0(*m_tmp_p0);

  result.copy_from(m_tmp);
}


} // end of namespace ocean
} // end of namespace pism
