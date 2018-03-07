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
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace ocean {

PicoGeometry::PicoGeometry(IceGrid::ConstPtr grid)
  : Component(grid),
    m_continental_shelf(grid, "pico_ocean_contshelf_mask", WITHOUT_GHOSTS),
    m_boxes(grid, "pico_ocean_box_mask", WITHOUT_GHOSTS),
    m_ice_shelves(grid, "pico_shelf_mask", WITHOUT_GHOSTS),
    m_ice_rises(grid, "ice_rises", WITHOUT_GHOSTS),
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
  compute_ice_rises(cell_type, m_ice_rises);

  compute_continental_shelf_mask(bed_elevation, m_ice_rises,
                                 m_config->get_double("ocean.pico.continental_shelf_depth"),
                                 m_continental_shelf);

  compute_ice_shelf_mask(m_ice_rises, m_ice_shelves);
}

/*!
 * Re-label components in a mask processed by label_connected_components.
 *
 * The biggest one gets the value of 2, all the other ones 1, the background is set to
 * zero.
 */
void PicoGeometry::relabel_by_size(IceModelVec2Int &mask) {

  int max_index = mask.range().max;

  if (max_index < 1) {
    // No components labeled. Fill the mask with zeros and quit.
    mask.set(0.0);
    return;
  }

  std::vector<double> area(max_index + 1, 0.0);
  {

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        int index = mask(i, j);

        if (index > max_index or index < 0) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "invalid component index: %d", index);
        }

        if (index > 0) {
          // count areas of actual components, ignoring the background (index == 0)
          area[index] += 1.0;
        }

      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

    for (unsigned int k = 0; k < area.size(); ++k) {
      area[k] = GlobalSum(m_grid->com, area[k]);
    }
  }

  // find the biggest component
  int biggest_component = 0;
  for (unsigned int k = 0; k < area.size(); ++k) {
    if (area[k] > area[biggest_component]) {
      biggest_component = k;
    }
  }

  // re-label
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int component_index = mask.as_int(i, j);

    if (component_index == biggest_component) {
      mask(i, j) = 2.0;
    } else if (component_index > 0) {
      mask(i, j) = 1.0;
    } else {
      mask(i, j) = 0.0;
    }
  }
}

/*!
 * Run the serial connected-component labeling algorithm on m_tmp.
 */
void PicoGeometry::label_tmp() {
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
}

/*!
 * Compute the mask identifying "subglacial lakes", i.e. floating ice areas that are not
 * connected to the open ocean.
 *
 * Resulting mask contains:
 *
 * 0 - grounded ice
 * 1 - floating ice not connected to the open ocean
 * 2 - floating ice or ice-free ocean connected to the open ocean
 */
void PicoGeometry::compute_lakes(const IceModelVec2CellType &cell_type,
                                 IceModelVec2Int &result) {
  IceModelVec::AccessList list{&cell_type, &m_tmp};

  // mask of zeros and ones: one if floating ice, zero otherwise
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j)) {
      m_tmp(i, j) = 1.0;
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  label_tmp();

  relabel_by_size(m_tmp);

  result.copy_from(m_tmp);
}

/*!
 * Compute the mask identifying "ice rises", i.e. grounded ice areas not connected to the
 * continental ice sheet.
 *
 * Resulting mask contains:
 *
 * 0 - ocean
 * 1 - ice rises
 * 2 - continental ice sheet
 * 3 - floating ice
 */
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

  label_tmp();

  relabel_by_size(m_tmp);

  // mark floating ice areas in this mask (reduces the number of masks we need later)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_tmp(i, j) == 0.0 and cell_type.icy(i, j)) {
      m_tmp(i, j) = FLOATING;
    }
  }

  result.copy_from(m_tmp);
}

/*!
 * Compute the continental ice shelf mask.
 *
 * Resulting mask contains:
 *
 * 0 - ocean or icy
 * 1 - ice-free areas with bed elevation > threshold and not connected to the continental ice sheet
 * 2 - ice-free areas with bed elevation > threshold, connected to the continental ice sheet
 */
void PicoGeometry::compute_continental_shelf_mask(const IceModelVec2S &bed_elevation,
                                                  const IceModelVec2Int &ice_rises_mask,
                                                  double bed_elevation_threshold,
                                                  IceModelVec2Int &result) {
  IceModelVec::AccessList list{&bed_elevation, &ice_rises_mask, &m_tmp};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_tmp(i, j) = 0.0;

    if (bed_elevation(i, j) > bed_elevation_threshold) {
      m_tmp(i, j) = 1.0;
    }

    if (ice_rises_mask.as_int(i, j) == CONTINENTAL) {
      m_tmp(i, j) = 2.0;
    }
  }

  // use "iceberg identification" to label parts *not* connected to the continental ice
  // sheet
  {
    m_tmp.put_on_proc0(*m_tmp_p0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        petsc::VecArray mask_p0(*m_tmp_p0);
        label_connected_components(mask_p0.get(), m_grid->My(), m_grid->Mx(), true, 2.0);
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();

    m_tmp.get_from_proc0(*m_tmp_p0);
  }

  // At this point areas with bed > threshold are 1, everything else is zero.
  //
  // Now we need to mark the continental shelf itself.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_tmp(i, j) > 0.0) {
      continue;
    }

    if (bed_elevation(i, j) > bed_elevation_threshold and
        ice_rises_mask.as_int(i, j) == OCEAN) {
      m_tmp(i, j) = 2.0;
    }
  }

  result.copy_from(m_tmp);
}

/*!
 * Compute the mask identifying ice shelves.
 *
 * Each shelf gets an individual integer label.
 *
 * Two shelves connected by an ice rise are considered to be parts of the same shelf.
 */
void PicoGeometry::compute_ice_shelf_mask(const IceModelVec2Int &ice_rises_mask,
                                          IceModelVec2Int &result) {
  IceModelVec::AccessList list{&ice_rises_mask, &m_tmp};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int M = ice_rises_mask.as_int(i, j);

    if (M == RISE or M == FLOATING) {
      m_tmp(i, j) = 1.0;
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  label_tmp();

  // remove ice rises
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_rises_mask.as_int(i, j) == RISE) {
      m_tmp(i, j) = 0.0;
    }
  }

  result.copy_from(m_tmp);
}

/*!
 * Compute the mask identifying ice-free ocean and "holes" in ice shelves.
 *
 * Resulting mask contains:
 *
 * - 0 - icy cells
 * - 1 - ice-free ocean which is not connected to the open ocean
 * - 2 - open ocean
 *
 */
void PicoGeometry::compute_ocean_mask(const IceModelVec2CellType &cell_type,
                                      IceModelVec2Int &result) {
  IceModelVec::AccessList list{&cell_type, &m_tmp};

  // mask of zeros and ones: one if ice-free ocean, zero otherwise
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free_ocean(i, j)) {
      m_tmp(i, j) = 1.0;
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  label_tmp();

  relabel_by_size(m_tmp);

  result.copy_from(m_tmp);
}

void PicoGeometry::compute_distances_gl(const IceModelVec2CellType &mask,
                                        const IceModelVec2Int &ocean_mask,
                                        const IceModelVec2Int &ice_rises,
                                        IceModelVec2Int &result) {

 // to find DistGL, 1 if floating and directly adjacent to a grounded cell
  double current_label = 1;

  IceModelVec::AccessList list{ &mask, &ice_rises, &ocean_mask, &result };

  bool exclude_rises = true;

  result.set(0);

  const int EXCLUDE = 1;
  const int INNER = 2;

  // Find the grounding line and the ice front and
  // set DistGL to 1 if ice shelf cell is next to the grounding line,
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    bool condition;
    if (exclude_rises) {
      condition =
          (mask(i, j) == MASK_FLOATING || ice_rises(i, j) == EXCLUDE || ocean_mask(i, j) == EXCLUDE);
    } else {
      condition = (mask(i, j) == MASK_FLOATING || ocean_mask(i, j) == EXCLUDE);
    }

    if (condition) { //if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

      // label the shelf cells adjacent to the grounding line with DistGL = 1
      bool neighbor_to_land;
      if (exclude_rises) {
        neighbor_to_land =
            (ice_rises(i, j + 1) == INNER || ice_rises(i, j - 1) == INNER ||
             ice_rises(i + 1, j) == INNER || ice_rises(i - 1, j) == INNER ||
             ice_rises(i + 1, j + 1) == INNER || ice_rises(i + 1, j - 1) == INNER ||
             ice_rises(i - 1, j + 1) == INNER || ice_rises(i - 1, j - 1) == INNER);
      } else {
        neighbor_to_land =
            (mask(i, j + 1) < MASK_FLOATING || mask(i, j - 1) < MASK_FLOATING || mask(i + 1, j) < MASK_FLOATING ||
             mask(i - 1, j) < MASK_FLOATING || mask(i + 1, j + 1) < MASK_FLOATING || mask(i + 1, j - 1) < MASK_FLOATING ||
             mask(i - 1, j + 1) < MASK_FLOATING || mask(i - 1, j - 1) < MASK_FLOATING);
      }

      if (neighbor_to_land) {
        // i.e. there is a grounded neighboring cell (which is not ice rise!)
        result(i, j) = current_label;
      } // no else

    }
  }

  result.update_ghosts();

  // DistGL calculation: Derive the distance from the grounding line for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  double continue_loop = 1;
  while (continue_loop != 0) {

    continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask(i, j) == MASK_FLOATING or
          ocean_mask(i, j) == EXCLUDE or
          (exclude_rises and ice_rises(i, j) == EXCLUDE)) {
        // this cell is floating or an hole in the ice shelf (or an ice rise)

        auto R = result.int_star(i, j);

        if (R.ij == 0 and
            (R.n == current_label or R.s == current_label or
             R.e == current_label or R.w == current_label)) {
          // i.e. this is an shelf cell with no distance assigned yet and with a neighbor
          // that has a distance assigned
          result(i, j) = current_label + 1;
          continue_loop = 1;
        }
      }
    } // loop over grid points

    current_label++;
    result.update_ghosts();

    continue_loop = GlobalMax(m_grid->com, continue_loop);
  }
}

void PicoGeometry::compute_distances_if(const IceModelVec2CellType &mask,
                                        const IceModelVec2Int &ocean_mask,
                                        const IceModelVec2Int &ice_rises,
                                        IceModelVec2Int &dist_if) {

  double currentLabelIF = 1; // to find DistIF, 1 if floating and directly adjacent to an ocean cell

  double global_continue_loop = 1;
  double local_continue_loop  = 0;

  IceModelVec::AccessList list{ &mask, &ice_rises, &ocean_mask, &dist_if };

  bool exclude_rises = true;

  dist_if.set(0);

  const int EXCLUDE = 1;
  const int INNER = 2;

  // Find the grounding line and the ice front and
  // set DistIF to 1 if ice shelf cell is next to the calving front.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    bool condition;
    if (exclude_rises) {
      condition =
          (mask(i, j) == MASK_FLOATING || ice_rises(i, j) == EXCLUDE || ocean_mask(i, j) == EXCLUDE);
    } else {
      condition = (mask(i, j) == MASK_FLOATING || ocean_mask(i, j) == EXCLUDE);
    }

    if (condition) { //if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

      // label the shelf cells adjacent to the calving front with DistIF = 1,
      // we do not need to exclude ice rises in this case.
      bool neighbor_to_ocean;
      neighbor_to_ocean = (ocean_mask(i, j + 1) == INNER || ocean_mask(i, j - 1) == INNER ||
                           ocean_mask(i + 1, j) == INNER || ocean_mask(i - 1, j) == INNER);

      if (neighbor_to_ocean) {
        dist_if(i, j) = currentLabelIF;
      }
    }
  }

  dist_if.update_ghosts();

  // DistIF calculation: Derive the distance from the calving front for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1; // start loop
  while (global_continue_loop != 0) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exclude_rises) {
        condition = (mask(i, j) == MASK_FLOATING || ice_rises(i, j) == EXCLUDE ||
                     ocean_mask(i, j) == EXCLUDE);
      } else {
        condition = (mask(i, j) == MASK_FLOATING || ocean_mask(i, j) == EXCLUDE);
      }

      if (condition && dist_if(i, j) == 0 &&
          (dist_if(i, j + 1) == currentLabelIF || dist_if(i, j - 1) == currentLabelIF ||
           dist_if(i + 1, j) == currentLabelIF || dist_if(i - 1, j) == currentLabelIF)) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
        dist_if(i, j) = currentLabelIF + 1;
        local_continue_loop = 1;
      } //if

    } // for

    currentLabelIF++;
    dist_if.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistIF
}
} // end of namespace ocean
} // end of namespace pism
