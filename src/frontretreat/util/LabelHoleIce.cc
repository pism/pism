/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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

#include "LabelHoleIce.hh"

#include "pism/util/connected_components.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace calving {

LabelHoleIce::LabelHoleIce(IceGrid::ConstPtr g)
  : Component(g),
    m_forced_open_ocean_mask(m_grid, "forced_open_ocean", WITHOUT_GHOSTS),
    m_enclosed_ocean_mask(m_grid, "enclosed_ocean_mask", WITHOUT_GHOSTS){

  m_mask_enclose_ocean_p0 = m_enclosed_ocean_mask.allocate_proc0_copy();
}

LabelHoleIce::~LabelHoleIce() {
  // empty
}

void LabelHoleIce::init() {
}


/**
 * Define points that are considered as open ocean in contrast to an enclosed ocean, which are holes in ice shelves.
 *
 * @param[in,out] forced open ocean masked
 */
void LabelHoleIce::open_ocean_mask1(const IceModelVec2T &retreat_mask,
				    const IceModelVec2S &bed,
				    const IceModelVec2S &sea_level,
				    IceModelVec2Int &open_ocean_mask) {
  const double depth_abyssal = 2000.0;
  const double depth_coast = 0.1;

  {
    IceModelVec::AccessList list{&retreat_mask, &bed, &sea_level, &open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      double depth_ocean = bed(i, j)+sea_level(i, j) //todo:??:;(*bed)(i, j)+(*sea_level)(i, j)???

	if (grid_edge(*m_grid, i, j) && depth_ocean < depth_abyssal) {
	   // Abyssal ocean at the domain margin
	   open_ocean_mask(i, j) = 1;
	} else if (retreat_mask(i, j) > 0.5 && depth_ocean < depth_coast) {
	   // Forced retreat
	   open_ocean_mask(i, j) = 1;
	} else {
	   // Otherwise
	   open_ocean_mask(i, j) = 0;
	}
    }
  }
}

/**
 * Define points that are considered as open ocean in contrast to an enclosed ocean, which are holes in ice shelves.
 *
 * @param[in,out] forced open ocean masked
 */
void LabelHoleIce::open_ocean_mask2(const IceModelVec2S &bed,
				    const IceModelVec2S &sea_level,
				    IceModelVec2Int &open_ocean_mask) {
  const double depth_abyssal = 2000.0;

  {
    IceModelVec::AccessList list{&bed, &sea_level, &open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      double depth_ocean = bed(i, j)+sea_level(i, j)

	if (grid_edge(*m_grid, i, j) && depth_ocean < depth_abyssal) {
	   // Abyssal ocean at the domain margin
	   open_ocean_mask(i, j) = 1;
	} else {
	   // Otherwise
	   open_ocean_mask(i, j) = 0;
	}
    }
  }
}


/**
 * Use PISM's open ocean mask to identify holes in ice shelves avoiding "black-hole" calving.
 *
 * @param[in,out] pism_mask PISM's ice cover mask
 */
void LabelHoleIce::update(const IceModelVec2Int &open_ocean_mask,
			  IceModelVec2CellType &mask) {
  const int
    mask_not_enclosed_ocean = 1,
    mask_enclosed_ocean = 2;

  // prepare the mask that will be handed to the connected component
  // labeling code:
  {
    m_enclosed_ocean_mask.set(0.0);

    IceModelVec::AccessList list{&open_ocean_maskm, &mask, &m_enclosed_ocean_mask};

    // Ocean points are potentially enclosed ocean points
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.ocean(i, j) == true) {
        m_enclosed_ocean_mask(i, j) = mask_enclosed_ocean;
      }
    }

    // Open ocean points are not enclosed and defined by open_ocean_mask=1.
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (open_ocean_mask(i, j) ==  1) {
        m_enclosed_ocean_mask(i, j) = mask_not_enclosed_ocean;
      }
    }
  }

  // identify holes in ice shelves using serial code on processor 0:
  {
    m_enclosed_ocean_mask.put_on_proc0(*m_mask_enclose_ocean_p0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        petsc::VecArray mask_p0(*m_mask_enclose_ocean_p0);
        label_connected_components(mask_p0.get(), m_grid->My(), m_grid->Mx(), true, mask_not_enclosed_ocean);
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();

    m_enclosed_ocean_mask.get_from_proc0(*m_mask_enclose_ocean_p0);
  }

  // correct ice thickness and the cell type mask using the resulting
  // "ice shelf hole" mask:
  {
    IceModelVec::AccessList list{&mask, &m_enclosed_ocean_mask, &open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_enclosed_ocean_mask(i,j) > 0.5 && open_ocean_mask(i,j) < 1) {
        mask(i,j)     = MASK_ICE_FREE_ENCLOSED_OCEAN;
      }
    }
  }

  // update ghosts of the mask
  mask.update_ghosts();
}

} // end of namespace calving
} // end of namespace pism
