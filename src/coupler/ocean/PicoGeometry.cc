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

#include <algorithm> // max_element

#include "PicoGeometry.hh"
#include "pism/calving/connected_components.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace ocean {

PicoGeometry::PicoGeometry(IceGrid::ConstPtr grid)
    : Component(grid),
      m_continental_shelf(grid, "pico_contshelf_mask", WITHOUT_GHOSTS),
      m_boxes(grid, "pico_box_mask", WITHOUT_GHOSTS),
      m_ice_shelves(grid, "pico_shelf_mask", WITHOUT_GHOSTS),
      m_distance_gl(grid, "pico_distance_gl", WITH_GHOSTS),
      m_distance_cf(grid, "pico_distance_cf", WITH_GHOSTS),
      m_ocean_mask(grid, "pico_ocean_mask", WITH_GHOSTS),
      m_lake_mask(grid, "pico_lake_mask", WITHOUT_GHOSTS),
      m_ice_rises(grid, "pico_ice_rise_mask", WITH_GHOSTS),
      m_tmp(grid, "temporary_storage", WITHOUT_GHOSTS) {

  m_boxes.metadata().set_double("_FillValue", 0.0);

  m_ice_rises.metadata().set_doubles("flag_values",
                                     {OCEAN, RISE, CONTINENTAL, FLOATING});
  m_ice_rises.metadata().set_string("flag_meanings",
                                     "ocean ice_rise continental_ice_sheet, floating_ice");

  m_tmp_p0 = m_tmp.allocate_proc0_copy();
}

PicoGeometry::~PicoGeometry() {
  // empty
}

const IceModelVec2Int &PicoGeometry::continental_shelf_mask() const {
  return m_continental_shelf;
}

const IceModelVec2Int &PicoGeometry::box_mask() const {
  return m_boxes;
}

const IceModelVec2Int &PicoGeometry::ice_shelf_mask() const {
  return m_ice_shelves;
}

const IceModelVec2Int &PicoGeometry::ice_rise_mask() const {
  return m_ice_rises;
}

/*!
 * Compute masks needed by the PICO physics code.
 *
 * After this call box_mask(), ice_shelf_mask(), and continental_shelf_mask() will be up
 * to date.
 */
void PicoGeometry::update(const IceModelVec2S &bed_elevation, const IceModelVec2CellType &cell_type) {
  bool exclude_ice_rises = m_config->get_boolean("ocean.pico.exclude_ice_rises");

  int n_boxes = m_config->get_double("ocean.pico.number_of_boxes");

  double continental_shelf_depth = m_config->get_double("ocean.pico.continental_shelf_depth");

  // these three could be done at the same time
  {
    compute_ice_rises(cell_type, exclude_ice_rises, m_ice_rises);

    compute_ocean_mask(cell_type, m_ocean_mask);

    compute_lakes(cell_type, m_lake_mask);
  }

  // these two could be optimized by trying to reduce the number of times we update ghosts
  {
    m_ice_rises.update_ghosts();
    m_ocean_mask.update_ghosts();

    compute_distances_gl(m_ocean_mask, m_ice_rises, exclude_ice_rises, m_distance_gl);

    compute_distances_cf(m_ocean_mask, m_ice_rises, exclude_ice_rises, m_distance_cf);
  }

  // these two could be done at the same time
  {
    compute_ice_shelf_mask(m_ice_rises, m_lake_mask, m_ice_shelves);

    compute_continental_shelf_mask(bed_elevation, m_ice_rises, continental_shelf_depth, m_continental_shelf);
  }

  compute_box_mask(m_distance_gl, m_distance_cf, m_ice_shelves, n_boxes, m_boxes);
}


enum RelabelingType {BY_AREA, AREA_THRESHOLD};

/*!
 * Re-label components in a mask processed by label_connected_components.
 *
 * If type is `BY_AREA`, the biggest one gets the value of 2, all the other ones 1, the
 * background is set to zero.
 *
 * If type is `AREA_THRESHOLD`, patches with areas above `threshold` get the value of 2,
 * all the other ones 1, the background is set to zero.
 */
static void relabel(RelabelingType type,
                    double threshold,
                    IceModelVec2Int &mask) {

  IceGrid::ConstPtr grid = mask.grid();

  int max_index = mask.range().max;

  if (max_index < 1) {
    // No components labeled. Fill the mask with zeros and quit.
    mask.set(0.0);
    return;
  }

  std::vector<double> area(max_index + 1, 0.0);
  {

    ParallelSection loop(grid->com);
    try {
      for (Points p(*grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        int index = mask(i, j);

        if (index > max_index or index < 0) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid component index: %d", index);
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
      area[k] = grid->cell_area() * GlobalSum(grid->com, area[k]);
    }
  }

  if (type == BY_AREA) {
    // find the biggest component
    int biggest_component = 0;
    for (unsigned int k = 0; k < area.size(); ++k) {
      if (area[k] > area[biggest_component]) {
        biggest_component = k;
      }
    }

    // re-label
    for (Points p(*grid); p; p.next()) {
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
  } else {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int component_index = mask.as_int(i, j);

      if (area[component_index] > threshold) {
        mask(i, j) = 2.0;
      } else if (component_index > 0) {
        mask(i, j) = 1.0;
      } else {
        mask(i, j) = 0.0;
      }
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

static bool edge_p(int i, int j, int Mx, int My) {
  return (i == 0) or (i == Mx - 1) or (j == 0) or (j == My - 1);
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
void PicoGeometry::compute_lakes(const IceModelVec2CellType &cell_type, IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &cell_type, &m_tmp };

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My();

  // assume that ocean points (i.e. floating, either icy or ice-free) at the edge of the
  // domain belong to the "open ocean"
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j)) {
      m_tmp(i, j) = 1.0;

      if (edge_p(i, j, Mx, My)) {
        m_tmp(i, j) = 2.0;
      }
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  // identify "floating" areas that are not connected to the open ocean as defined above
  {
    m_tmp.put_on_proc0(*m_tmp_p0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        petsc::VecArray mask_p0(*m_tmp_p0);
        label_connected_components(mask_p0.get(), My, Mx, true, 2.0);
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();

    m_tmp.get_from_proc0(*m_tmp_p0);
  }

  result.copy_from(m_tmp);
}

/*!
 * Compute the mask identifying ice rises, i.e. grounded ice areas not connected to the
 * continental ice sheet.
 *
 * Resulting mask contains:
 *
 * 0 - ocean
 * 1 - ice rises
 * 2 - continental ice sheet
 * 3 - floating ice
 */
void PicoGeometry::compute_ice_rises(const IceModelVec2CellType &cell_type, bool exclude_ice_rises,
                                     IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &cell_type, &m_tmp };

  // mask of zeros and ones: one if grounded ice, zero otherwise
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.grounded(i, j)) {
      m_tmp(i, j) = 2.0;
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  if (exclude_ice_rises) {
    label_tmp();

    relabel(AREA_THRESHOLD,
            m_config->get_double("ocean.pico.maximum_ice_rise_area", "m2"),
            m_tmp);
  }

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
                                                  const IceModelVec2Int &ice_rise_mask, double bed_elevation_threshold,
                                                  IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &bed_elevation, &ice_rise_mask, &m_tmp };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_tmp(i, j) = 0.0;

    if (bed_elevation(i, j) > bed_elevation_threshold) {
      m_tmp(i, j) = 1.0;
    }

    if (ice_rise_mask.as_int(i, j) == CONTINENTAL) {
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
        ice_rise_mask.as_int(i, j) == OCEAN) {
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
 *
 * Floating ice cells that are not connected to the ocean ("subglacial lakes") are
 * excluded.
 */
void PicoGeometry::compute_ice_shelf_mask(const IceModelVec2Int &ice_rise_mask, const IceModelVec2Int &lake_mask,
                                          IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &ice_rise_mask, &lake_mask, &m_tmp };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int M = ice_rise_mask.as_int(i, j);

    if (M == RISE or M == FLOATING) {
      m_tmp(i, j) = 1.0;
    } else {
      m_tmp(i, j) = 0.0;
    }
  }

  label_tmp();

  // remove ice rises and lakes
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_rise_mask.as_int(i, j) == RISE or lake_mask.as_int(i, j) == 1) {
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
void PicoGeometry::compute_ocean_mask(const IceModelVec2CellType &cell_type, IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &cell_type, &m_tmp };

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

  relabel(BY_AREA, 0.0, m_tmp);

  result.copy_from(m_tmp);
}

void PicoGeometry::compute_distances_gl(const IceModelVec2Int &ocean_mask,
                                        const IceModelVec2Int &ice_rises,
                                        bool exclude_ice_rises,
                                        IceModelVec2Int &result) {

  IceModelVec::AccessList list{ &ice_rises, &ocean_mask, &result };

  result.set(-1);

  // Find the grounding line and the ice front and set result to 1 if ice shelf cell is
  // next to the grounding line, Ice holes within the shelf are treated like ice shelf
  // cells, if exclude_ice_rises is set then ice rises are also treated as ice shelf cells.

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (ice_rises.as_int(i, j) == FLOATING or
          ocean_mask.as_int(i, j) == 1 or
          (exclude_ice_rises and ice_rises.as_int(i, j) == RISE)) {
        // if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

        // label the shelf cells adjacent to the grounding line with 1
        bool neighbor_to_land =
          (ice_rises(i, j + 1) == CONTINENTAL or
           ice_rises(i, j - 1) == CONTINENTAL or
           ice_rises(i + 1, j) == CONTINENTAL or
           ice_rises(i - 1, j) == CONTINENTAL or
           ice_rises(i + 1, j + 1) == CONTINENTAL or
           ice_rises(i + 1, j - 1) == CONTINENTAL or
           ice_rises(i - 1, j + 1) == CONTINENTAL or
           ice_rises(i - 1, j - 1) == CONTINENTAL);

        if (neighbor_to_land) {
          // i.e. there is a grounded neighboring cell (which is not an ice rise!)
          result(i, j) = 1;
        } else {
          result(i, j) = 0;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();

  eikonal_equation(result);
}

void PicoGeometry::compute_distances_cf(const IceModelVec2Int &ocean_mask,
                                        const IceModelVec2Int &ice_rises,
                                        bool exclude_ice_rises,
                                        IceModelVec2Int &result) {

  IceModelVec::AccessList list{ &ice_rises, &ocean_mask, &result };

  result.set(-1);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (ice_rises.as_int(i, j) == FLOATING or
          ocean_mask.as_int(i, j) == 1 or
          (exclude_ice_rises and ice_rises.as_int(i, j) == RISE)) {
        // if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

        // label the shelf cells adjacent to the ice front with 1
        auto M = ocean_mask.int_star(i, j);

        if (M.n == 2 or M.e == 2 or M.s == 2 or M.w == 2) {
          // i.e. there is a neighboring open ocean cell
          result(i, j) = 1;
        } else {
          result(i, j) = 0;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();

  eikonal_equation(result);
}

/*!
 * Find an approximate solution of the Eikonal equation on a given domain.
 *
 * To specify the problem, the input field (mask) should be filled with
 *
 * - negative values outside the domain
 * - zeros within the domain
 * - ones at "wave front" locations
 *
 * For example, to compute distances from the grounding line within ice shelves, fill
 * generic ice shelf locations with zeros, set neighbors of the grounding line to 1, and
 * the rest of the grid with -1 or some other negative number.
 *
 * Note that this implementation updates ghosts *every* iteration. We could speed this
 * up by checking if a point at a boundary of the processor sub-domain was updated and
 * update ghosts in those cases only.
 */
void eikonal_equation(IceModelVec2Int &mask) {

  assert(mask.stencil_width() > 0);

  IceGrid::ConstPtr grid = mask.grid();

  double current_label = 1;
  double continue_loop = 1;
  while (continue_loop != 0) {

    continue_loop = 0;

    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.as_int(i, j) == 0) {

        auto R = mask.int_star(i, j);

        if (R.ij == 0 and
            (R.n == current_label or R.s == current_label or
             R.e == current_label or R.w == current_label)) {
          // i.e. this is an shelf cell with no distance assigned yet and with a neighbor
          // that has a distance assigned
          mask(i, j) = current_label + 1;
          continue_loop = 1;
        }
      }
    } // loop over grid points

    current_label++;
    mask.update_ghosts();

    continue_loop = GlobalMax(grid->com, continue_loop);
  }
}

void PicoGeometry::compute_box_mask(const IceModelVec2Int &D_gl, const IceModelVec2Int &D_cf,
                                    const IceModelVec2Int &shelf_mask, int max_number_of_boxes,
                                    IceModelVec2Int &result) {

  IceModelVec::AccessList list{ &D_gl, &D_cf, &shelf_mask, &result };

  int n_shelves = shelf_mask.range().max + 1;

  std::vector<double> GL_distance_max(n_shelves, 0.0);
  std::vector<double> CF_distance_max(n_shelves, 0.0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask(i, j);
    assert(shelf_id >= 0);
    assert(shelf_id < n_shelves + 1);

    if (shelf_id == 0) {
      // not at a shelf; skip to the next grid point
      continue;
    }

    if (D_gl(i, j) > GL_distance_max[shelf_id]) {
      GL_distance_max[shelf_id] = D_gl(i, j);
    }

    if (D_cf(i, j) > CF_distance_max[shelf_id]) {
      CF_distance_max[shelf_id] = D_cf(i, j);
    }
  }

  // compute global maximums
  for (int k = 0; k < n_shelves; ++k) {
    GL_distance_max[k] = GlobalMax(m_grid->com, GL_distance_max[k]);
    CF_distance_max[k] = GlobalMax(m_grid->com, CF_distance_max[k]);
  }

  double GL_distance_ref = *std::max_element(GL_distance_max.begin(), GL_distance_max.end());

  // compute the number of boxes in each shelf

  std::vector<int> n_boxes(n_shelves, 0);
  int n_min   = 1;
  double zeta = 0.5;

  for (int k = 0; k < n_shelves; ++k) {
    n_boxes[k] = n_min + round(pow((GL_distance_max[k] / GL_distance_ref), zeta) * (max_number_of_boxes - n_min));

    n_boxes[k] = std::min(n_boxes[k], max_number_of_boxes);
  }

  result.set(0.0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int d_gl = D_gl.as_int(i, j);
    int d_cf = D_cf.as_int(i, j);

    if (shelf_mask.as_int(i, j) > 0 and d_gl > 0.0 and d_cf > 0.0) {
      int shelf_id = shelf_mask(i, j);
      int n = n_boxes[shelf_id];

      // relative position on the shelf (ranges from 0 to 1), increasing towards the
      // calving front (see equation 10 in the PICO paper)
      double r = d_gl / (double)(d_gl + d_cf);

      double C = pow(1.0 - r, 2.0);
      for (int k = 0; k < n; ++k) {
        if ((n - k - 1) / (double)n <= C and C <= (n - k) / (double)n) {
          result(i, j) = std::min(d_gl, k + 1);
        }
      }
    }
  }
}

} // end of namespace ocean
} // end of namespace pism
