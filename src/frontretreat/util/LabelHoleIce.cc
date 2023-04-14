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

#include "pism/frontretreat/util/LabelHoleIce.hh"
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
    m_bc_open_ocean_mask(m_grid, "bc_open_ocean_mask", WITHOUT_GHOSTS),
    m_enclosed_ocean_mask(m_grid, "enclosed_ocean_mask", WITHOUT_GHOSTS)
{
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
 */
void LabelHoleIce::open_ocean_mask_margin_retreat(const IceModelVec2S &bed,
						  const IceModelVec2S &sea_level,
						  const IceModelVec2S &ice_area_specific_volume,
						  const IceModelVec2S &ice_thickness) {
  //todo:alternative;const IceModelVec2S &retreat_mask instead of "ice_area_specific_volume/ice_thickness"
  const double depth_abyssal = 2000.0; //FIXME:changeable parameter
  const double depth_coast = 0.1;      //FIXME:changeable parameter

  m_bc_open_ocean_mask.set(0); //todo:erase?
  {
    IceModelVec::AccessList list{&bed, &sea_level, &ice_area_specific_volume, &ice_thickness, &m_bc_open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      double depth_ocean = bed(i, j)+sea_level(i, j);
      double retreat_factor;

      //todo:rm!;if ( ice_thickness(i, j) > 0.001 ) {
      //todo:rm!;	retreat_factor = ice_area_specific_volume(i, j)/ice_thickness(i, j);
      //todo:rm!;} else {
      //todo:rm!;	retreat_factor = std::min(ice_area_specific_volume(i, j), 1.0);
      //todo:rm!;}

      retreat_factor = ice_area_specific_volume(i, j)/std::max(0.0001, ice_thickness(i, j));
      //todo:alternative; retreat_factor = retreat_mask(i, j)

      if (grid_edge(*m_grid, i, j) && depth_ocean > depth_abyssal) {
	   // Abyssal ocean at the domain edge
	   m_bc_open_ocean_mask(i, j) = 1;
	} else if (retreat_factor > 0.5 && depth_ocean > depth_coast) {
	   // Forced retreat
	   m_bc_open_ocean_mask(i, j) = 1;
	} else {
	   // Otherwise
	   m_bc_open_ocean_mask(i, j) = 0;
	}
    }
  }
}

/**
 * Define points that are considered as open ocean in contrast to an enclosed ocean, which are holes in ice shelves.
 *
 */
void LabelHoleIce::open_ocean_mask_margin(const IceModelVec2S &bed,
					  const IceModelVec2S &sea_level) {

  const double depth_abyssal = 2000.0; //FIXME:changeable parameter

  m_bc_open_ocean_mask.set(0); //todo:erase?
  {
    IceModelVec::AccessList list{&bed, &sea_level, &m_bc_open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      double depth_ocean = bed(i, j)+sea_level(i, j);

	if (grid_edge(*m_grid, i, j) && depth_ocean > depth_abyssal) {
	   // Abyssal ocean at the domain margin
	   m_bc_open_ocean_mask(i, j) = 1;
	} else {
	   // Otherwise
	   m_bc_open_ocean_mask(i, j) = 0;
	}
    }
  }
}

/**
 * Use PISM's open ocean mask to identify holes in ice shelves avoiding "black-hole" calving.
 *
 * @param[in,out] pism_mask PISM's ice cover mask
 */
void LabelHoleIce::update(IceModelVec2CellType &mask) {
  const int
    mask_not_enclosed_ocean = 1,
    mask_enclosed_ocean = 2;

  // prepare the mask that will be handed to the connected component
  // labeling code:
  {
    m_enclosed_ocean_mask.set(0.0);

    IceModelVec::AccessList list{&m_bc_open_ocean_mask, &mask, &m_enclosed_ocean_mask};

    // Ocean points are potentially enclosed ocean points
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.ocean(i, j) == true) {
        m_enclosed_ocean_mask(i, j) = mask_enclosed_ocean;
      }
    }

    // Open ocean points are not enclosed and defined by m_bc_open_ocean_mask=1.
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_bc_open_ocean_mask(i, j) > 0.5 and mask.ocean(i, j)) { //todo:rm;if (m_bc_open_ocean_mask(i, j) ==  1) {
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
    IceModelVec::AccessList list{&mask, &m_enclosed_ocean_mask, &m_bc_open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_enclosed_ocean_mask(i, j) > 0.5 && m_bc_open_ocean_mask(i, j) < 0.5) {
        mask(i, j) = MASK_ICE_FREE_ENCLOSED_OCEAN;
      }
    }
  }

  // update ghosts of the mask
  mask.update_ghosts();
}

} // end of namespace calving
} // end of namespace pism
