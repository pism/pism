/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2023 PISM Authors
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
#include "pism/util/Diagnostic.hh"

namespace pism {
namespace calving {

LabelHoleIce::LabelHoleIce(IceGrid::ConstPtr g)
  : Component(g),
    m_bc_open_ocean_mask(m_grid, "bc_open_ocean_mask", WITHOUT_GHOSTS),
    m_enclosed_ocean_mask(m_grid, "enclosed_ocean_mask", WITHOUT_GHOSTS)
{
  m_bc_open_ocean_mask.set_attrs("diagnostic",
                                 "hole labelling seed mask: 1 where a cell is treated as"
                                 " open ocean a priori (abyssal ocean at the domain edge,"
                                 " or a forced-retreat point)",
                                 "", "", "", 0);
  m_enclosed_ocean_mask.set_attrs("diagnostic",
                                  "hole mask: 1 where an ocean cell was identified as"
                                  " enclosed by ice or land (a hole in an ice shelf),"
                                  " 0 otherwise",
                                  "", "", "", 0);

  // This object is allocated even when geometry.label_holes is false, so that
  // its diagnostics can always be named in -extra_vars. In that case update()
  // is never called and these fields are never assigned, so zero them here:
  // otherwise the reported values would be whatever happened to be in memory.
  m_bc_open_ocean_mask.set(0.0);
  m_enclosed_ocean_mask.set(0.0);

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

  {
    IceModelVec::AccessList list{&bed, &sea_level, &ice_area_specific_volume, &ice_thickness, &m_bc_open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      // Water depth, positive downward. "bed" is an elevation (negative below
      // sea level), so the depth_abyssal/depth_coast thresholds below only make
      // sense this way round; adding the two gives a quantity that is negative
      // everywhere in the ocean and never exceeds depth_abyssal.
      double depth_ocean = sea_level(i, j) - bed(i, j);
      double retreat_factor;

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

  {
    IceModelVec::AccessList list{&bed, &sea_level, &m_bc_open_ocean_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      // Water depth, positive downward. "bed" is an elevation (negative below
      // sea level), so the depth_abyssal/depth_coast thresholds below only make
      // sense this way round; adding the two gives a quantity that is negative
      // everywhere in the ocean and never exceeds depth_abyssal.
      double depth_ocean = sea_level(i, j) - bed(i, j);

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

    // Ice-free ocean points are potentially enclosed ocean points.
    //
    // NOTE: this must be ice_free_ocean(), not ocean(). ocean() is true for
    // floating ice too, which would make the ice shelf itself traversable by
    // the connected-component search below: a hole would then always be
    // "connected" to the open sea through the surrounding shelf and could never
    // be identified as enclosed. It would also write MASK_ICE_FREE_ENCLOSED_OCEAN
    // onto floating-ice cells, and mask 5 is ice-free as far as every predicate
    // is concerned.
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.ice_free_ocean(i, j) == true) {
        m_enclosed_ocean_mask(i, j) = mask_enclosed_ocean;
      }
    }

    // Open ocean points are not enclosed and defined by m_bc_open_ocean_mask=1.
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_bc_open_ocean_mask(i, j) > 0.5 and mask.ice_free_ocean(i, j)) {
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

/*!
 * Report both stages of the hole labelling. MASK_ICE_FREE_ENCLOSED_OCEAN itself
 * never reaches the output file, because Geometry::ensure_consistency()
 * recomputes cell_type after calving and GeometryCalculator only emits
 * 0/2/3/4. Reporting these two fields instead separates the seeding from the
 * connected-component labelling, so a failure in either stage is visible:
 *
 *   bc_open_ocean_mask  - which cells seeded the "this is open ocean" search
 *   enclosed_ocean_mask - which ocean cells came out as holes
 */
DiagnosticList LabelHoleIce::diagnostics_impl() const {
  return {{"bc_open_ocean_mask",  Diagnostic::wrap(m_bc_open_ocean_mask)},
          {"enclosed_ocean_mask", Diagnostic::wrap(m_enclosed_ocean_mask)}};
}

} // end of namespace calving
} // end of namespace pism
