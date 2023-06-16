/* Copyright (C) 2018, 2019, 2020, 2021, 2022, 2023 PISM Authors
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
#include "pism/util/connected_components.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/petscwrappers/Vec.hh"

#include "pism/coupler/util/options.hh"

namespace pism {
namespace ocean {

PicoGeometry::PicoGeometry(std::shared_ptr<const Grid> grid)
    : Component(grid),
      m_continental_shelf(grid, "pico_contshelf_mask"),
      m_boxes(grid, "pico_box_mask"),
      m_ice_shelves(grid, "pico_shelf_mask"),
      m_basin_mask(m_grid, "basins"),
      m_distance_gl(grid, "pico_distance_gl"),
      m_distance_cf(grid, "pico_distance_cf"),
      m_ocean_mask(grid, "pico_ocean_mask"),
      m_lake_mask(grid, "pico_lake_mask"),
      m_ice_rises(grid, "pico_ice_rise_mask"),
      m_tmp(grid, "temporary_storage") {

  m_continental_shelf.set_interpolation_type(NEAREST);
  m_boxes.set_interpolation_type(NEAREST);
  m_ice_shelves.set_interpolation_type(NEAREST);
  m_basin_mask.set_interpolation_type(NEAREST);
  m_distance_gl.set_interpolation_type(NEAREST);
  m_distance_cf.set_interpolation_type(NEAREST);
  m_ocean_mask.set_interpolation_type(NEAREST);
  m_lake_mask.set_interpolation_type(NEAREST);
  m_ice_rises.set_interpolation_type(NEAREST);
  m_tmp.set_interpolation_type(NEAREST);

  m_boxes.metadata()["_FillValue"] = {0.0};

  m_ice_rises.metadata()["flag_values"] = {OCEAN, RISE, CONTINENTAL, FLOATING};
  m_ice_rises.metadata()["flag_meanings"] =
    "ocean ice_rise continental_ice_sheet, floating_ice";

  m_basin_mask.set_attrs("climate_forcing", "mask determines basins for PICO",
                         "", "", "", 0);
  m_n_basins = 0;

  m_tmp_p0 = m_tmp.allocate_proc0_copy();
}

const array::Scalar &PicoGeometry::continental_shelf_mask() const {
  return m_continental_shelf;
}

const array::Scalar &PicoGeometry::box_mask() const {
  return m_boxes;
}

const array::Scalar &PicoGeometry::ice_shelf_mask() const {
  return m_ice_shelves;
}

const array::Scalar &PicoGeometry::ice_rise_mask() const {
  return m_ice_rises;
}

const array::Scalar &PicoGeometry::basin_mask() const {
  return m_basin_mask;
}

void PicoGeometry::init() {

  ForcingOptions opt(*m_grid->ctx(), "ocean.pico");

  m_basin_mask.regrid(opt.filename, CRITICAL);

  m_n_basins = static_cast<int>(max(m_basin_mask)) + 1;
}

/*!
 * Compute masks needed by the PICO physics code.
 *
 * After this call box_mask(), ice_shelf_mask(), and continental_shelf_mask() will be up
 * to date.
 */
void PicoGeometry::update(const array::Scalar &bed_elevation,
                          const array::CellType1 &cell_type) {

  // Update basin adjacency.
  //
  // basin_neighbors() below uses the cell type mask to find
  // adjacent basins by iterating over the current ice front. This means that basin
  // adjacency cannot be pre-computed during initialization.
  {
    m_basin_neighbors = basin_neighbors(cell_type, m_basin_mask);

    // report
    for (const auto &p : m_basin_neighbors) {
      std::vector<std::string> neighbors;
      for (const auto &n : p.second) {
        neighbors.emplace_back(pism::printf("%d", n));
      }
      std::string neighbor_list = pism::join(neighbors, ", ");
      m_log->message(3, "PICO: basin %d neighbors: %s\n",
                     p.first, neighbor_list.c_str());
    }
  }

  bool exclude_ice_rises = m_config->get_flag("ocean.pico.exclude_ice_rises");

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

  // computing ice_shelf_mask and box_mask could be done at the same time
  {
    compute_ice_shelf_mask(m_ice_rises, m_lake_mask, m_ice_shelves);
    auto n_shelves = static_cast<int>(max(m_ice_shelves)) + 1;

    std::vector<int> cfs_in_basins_per_shelf(n_shelves*m_n_basins, 0);
    std::vector<int> most_shelf_cells_in_basin(n_shelves, 0);
    identify_calving_front_connection(cell_type, m_basin_mask, m_ice_shelves, n_shelves,
                                      most_shelf_cells_in_basin, cfs_in_basins_per_shelf);

    split_ice_shelves(cell_type, m_basin_mask, m_basin_neighbors,
                      most_shelf_cells_in_basin, cfs_in_basins_per_shelf, n_shelves,
                      m_ice_shelves);

    double continental_shelf_depth = m_config->get_number("ocean.pico.continental_shelf_depth");

    compute_continental_shelf_mask(bed_elevation, m_ice_rises, continental_shelf_depth,
                                   m_continental_shelf);
  }

  int n_boxes = static_cast<int>(m_config->get_number("ocean.pico.number_of_boxes"));

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
                    array::Scalar &mask) {

  auto grid = mask.grid();

  int max_index = static_cast<int>(array::max(mask));

  if (max_index < 1) {
    // No components labeled. Fill the mask with zeros and quit.
    mask.set(0.0);
    return;
  }

  std::vector<double> area(max_index + 1, 0.0);
  std::vector<double> area1(max_index + 1, 0.0);
  {

    ParallelSection loop(grid->com);
    try {
      for (auto p = grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();

        int index = mask.as_int(i, j);

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
    GlobalSum(grid->com, area.data(), area1.data(), area.size());

    // copy data
    area = area1;

    for (unsigned int k = 0; k < area.size(); ++k) {
      area[k] = grid->cell_area() * area[k];
    }
  }

  if (type == BY_AREA) {
    // find the biggest component
    int biggest_component = 0;
    for (unsigned int k = 0; k < area.size(); ++k) {
      if (area[k] > area[biggest_component]) {
        biggest_component = static_cast<int>(k);
      }
    }

    // re-label
    for (auto p = grid->points(); p; p.next()) {
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
    for (auto p = grid->points(); p; p.next()) {
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
      label_connected_components(mask_p0.get(),
                                 static_cast<int>(m_grid->My()),
                                 static_cast<int>(m_grid->Mx()),
                                 false,
                                 0.0);
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
void PicoGeometry::compute_lakes(const array::CellType &cell_type, array::Scalar &result) {
  array::AccessScope list{ &cell_type, &m_tmp };

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My();

  // assume that ocean points (i.e. floating, either icy or ice-free) at the edge of the
  // domain belong to the "open ocean"
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j)) {
      m_tmp(i, j) = 1.0;

      if (grid::domain_edge(*m_grid, i, j)) {
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
void PicoGeometry::compute_ice_rises(const array::CellType &cell_type, bool exclude_ice_rises,
                                     array::Scalar &result) {
  array::AccessScope list{ &cell_type, &m_tmp };

  // mask of zeros and ones: one if grounded ice, zero otherwise
  for (auto p = m_grid->points(); p; p.next()) {
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
            m_config->get_number("ocean.pico.maximum_ice_rise_area", "m2"),
            m_tmp);
  }

  // mark floating ice areas in this mask (reduces the number of masks we need later)
  for (auto p = m_grid->points(); p; p.next()) {
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
void PicoGeometry::compute_continental_shelf_mask(const array::Scalar &bed_elevation,
                                                  const array::Scalar &ice_rise_mask,
                                                  double bed_elevation_threshold,
                                                  array::Scalar &result) {
  array::AccessScope list{ &bed_elevation, &ice_rise_mask, &m_tmp };

  for (auto p = m_grid->points(); p; p.next()) {
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
  for (auto p = m_grid->points(); p; p.next()) {
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
void PicoGeometry::compute_ice_shelf_mask(const array::Scalar &ice_rise_mask, const array::Scalar &lake_mask,
                                          array::Scalar &result) {
  array::AccessScope list{ &ice_rise_mask, &lake_mask, &m_tmp };

  for (auto p = m_grid->points(); p; p.next()) {
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
  for (auto p = m_grid->points(); p; p.next()) {
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
void PicoGeometry::compute_ocean_mask(const array::CellType &cell_type, array::Scalar &result) {
  array::AccessScope list{ &cell_type, &m_tmp };

  // mask of zeros and ones: one if ice-free ocean, zero otherwise
  for (auto p = m_grid->points(); p; p.next()) {
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

/*!
 * Find neighboring basins by checking for the basin boundaries on the ice free ocean.
 *
 * Should we identify the intersection at the coastline instead?
 *
 * Returns the map from the basin index to a set of indexes of neighbors.
 */
std::map<int,std::set<int> > PicoGeometry::basin_neighbors(const array::CellType1 &cell_type,
                                                           const array::Scalar1 &basin_mask) {
  using mask::ice_free_ocean;

  // Allocate the adjacency matrix. This uses twice the amount of storage necessary (the
  // matrix is symmetric), but in known cases (i.e. with around 20 basins) we're wasting
  // only ~200*sizeof(int) bytes, which is not that bad (and the code is a bit simpler).
  std::vector<int> adjacency_matrix(m_n_basins * m_n_basins, 0);

  // short-cuts
  auto mark_as_neighbors = [&](int b1, int b2) {
    adjacency_matrix[b1 * m_n_basins + b2] = 1;
    // preserve symmetry:
    adjacency_matrix[b2 * m_n_basins + b1] = 1;
  };

  auto adjacent = [&](int b1, int b2) {
    return adjacency_matrix[b1 * m_n_basins + b2] > 0;
  };

  array::AccessScope list{ &cell_type, &basin_mask };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto B = basin_mask.star(i, j);

    bool next_to_icefront = (cell_type.ice_free_ocean(i, j) and cell_type.next_to_ice(i,j));

    // skip the "dummy" basin and cells that are not at the ice front
    if (B.c == 0 or not next_to_icefront) {
      continue;
    }

    // Zero out IDs of basins for cell neighbors that are outside the modeling domain.
    //
    // This prevents "wrap-around" at grid boundaries.
    {
      B.n *= static_cast<int>(j < (int)m_grid->My() - 1);
      B.e *= static_cast<int>(i < (int)m_grid->Mx() - 1);
      B.s *= static_cast<int>(j > 0);
      B.w *= static_cast<int>(i > 0);
    }

    auto M = cell_type.star_int(i, j);

    if (ice_free_ocean(M.n)) {
      mark_as_neighbors(B.c, B.n);
    }

    if (ice_free_ocean(M.s)) {
      mark_as_neighbors(B.c, B.s);
    }

    if (ice_free_ocean(M.e)) {
      mark_as_neighbors(B.c, B.e);
    }

    if (ice_free_ocean(M.w)) {
      mark_as_neighbors(B.c, B.w);
    }
  }

  // Make a copy to allreduce efficiently with IntelMPI:
  {
    std::vector<int> tmp(adjacency_matrix.size(), 0);
    GlobalMax(m_grid->com, adjacency_matrix.data(), tmp.data(),
              static_cast<int>(tmp.size()));
    // Copy results:
    adjacency_matrix = tmp;
  }

  // Convert the matrix into a map "basin ID -> set of neighbors' IDs":
  std::map<int,std::set<int> > result;
  for (int b1 = 1; b1 < m_n_basins; ++b1) {
    for (int b2 = b1 + 1; b2 < m_n_basins; ++b2) {
      if (adjacent(b1, b2)) {
        result[b1].insert(b2);
        result[b2].insert(b1);
      }
    }
  }

  return result;
}

/*!
 *  Find all basins b, in which the ice shelf s has a calving front with potential ocean water intrusion.
 *  Find the basin bmax, in which the ice shelf s has the most cells
 */
void PicoGeometry::identify_calving_front_connection(const array::CellType1 &cell_type,
                                                     const array::Scalar &basin_mask,
                                                     const array::Scalar &shelf_mask,
                                                     int n_shelves,
                                                     std::vector<int> &most_shelf_cells_in_basin,
                                                     std::vector<int> &cfs_in_basins_per_shelf) {

  std::vector<int> n_shelf_cells_per_basin(n_shelves * m_n_basins,0);
  // additional vectors to allreduce efficiently with IntelMPI
  std::vector<int> n_shelf_cells_per_basinr(n_shelves * m_n_basins,0);
  std::vector<int> cfs_in_basins_per_shelfr(n_shelves * m_n_basins,0);
  std::vector<int> most_shelf_cells_in_basinr(n_shelves, 0);

  array::AccessScope list{ &cell_type, &basin_mask, &shelf_mask };

  {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();
      int s = shelf_mask.as_int(i, j);
      int b = basin_mask.as_int(i, j);
      int sb = s * m_n_basins + b;
      n_shelf_cells_per_basin[sb]++;

      if (cell_type.as_int(i, j) == MASK_FLOATING) {
        auto M = cell_type.star(i, j);
        if (M.n == MASK_ICE_FREE_OCEAN or
            M.e == MASK_ICE_FREE_OCEAN or
            M.s == MASK_ICE_FREE_OCEAN or
            M.w == MASK_ICE_FREE_OCEAN) {
          if (cfs_in_basins_per_shelf[sb] != b) {
            cfs_in_basins_per_shelf[sb] = b;
          }
        }
      }
    }

    GlobalSum(m_grid->com, cfs_in_basins_per_shelf.data(),
              cfs_in_basins_per_shelfr.data(), n_shelves*m_n_basins);
    GlobalSum(m_grid->com, n_shelf_cells_per_basin.data(),
              n_shelf_cells_per_basinr.data(), n_shelves*m_n_basins);
    // copy values
    cfs_in_basins_per_shelf = cfs_in_basins_per_shelfr;
    n_shelf_cells_per_basin = n_shelf_cells_per_basinr;

    for (int s = 0; s < n_shelves; s++) {
      int n_shelf_cells_per_basin_max = 0;
      for (int b = 0; b < m_n_basins; b++) {
        int sb = s * m_n_basins + b;
        if (n_shelf_cells_per_basin[sb] > n_shelf_cells_per_basin_max) {
          most_shelf_cells_in_basin[s] = b;
          n_shelf_cells_per_basin_max  = n_shelf_cells_per_basin[sb];
        }
      }
    }
  }
}


/*!
 * Find all ice shelves s that spread across non-neighboring basins with calving fronts in those basins and add an ice shelf mask number.
 */
void PicoGeometry::split_ice_shelves(const array::CellType &cell_type,
                                     const array::Scalar &basin_mask,
                                     const std::map<int, std::set<int> > &basin_neighbors,
                                     const std::vector<int> &most_shelf_cells_in_basin,
                                     const std::vector<int> &cfs_in_basins_per_shelf,
                                     int n_shelves,
                                     array::Scalar &shelf_mask) {
  m_tmp.copy_from(shelf_mask);

  std::vector<int> n_shelf_cells_to_split(n_shelves * m_n_basins, 0);

  array::AccessScope list{ &cell_type, &basin_mask, &shelf_mask, &m_tmp };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();
      if (cell_type.as_int(i, j) == MASK_FLOATING) {
        int b = basin_mask.as_int(i, j);
        int s = shelf_mask.as_int(i, j);
        int b0 = most_shelf_cells_in_basin[s];
        // basin_neighbors.at(b) may throw
        bool neighbors = basin_neighbors.at(b).count(b0) > 0;
        if (b != b0 and (not neighbors) and
            cfs_in_basins_per_shelf[s * m_n_basins + b] > 0) {
          n_shelf_cells_to_split[s * m_n_basins + b]++;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  {
    std::vector<int> tmp(n_shelves * m_n_basins, 0);
    GlobalSum(m_grid->com, n_shelf_cells_to_split.data(),
              tmp.data(), n_shelves * m_n_basins);
    // copy values
    n_shelf_cells_to_split = tmp;
  }

  // no GlobalSum needed here, only local:
  std::vector<int> add_shelf_instance(n_shelves * m_n_basins, 0);
  int n_shelf_numbers_to_add = 0;
  for (int s = 0; s < n_shelves; s++) {
    int b0 = most_shelf_cells_in_basin[s];
    for (int b = 0; b < m_n_basins; b++) {
      if (n_shelf_cells_to_split[s * m_n_basins + b] > 0) {
        n_shelf_numbers_to_add += 1;
        add_shelf_instance[s * m_n_basins + b] = n_shelves + n_shelf_numbers_to_add;
        m_log->message(3, "\nPICO, split ice shelf s=%d with bmax=%d "
                       "and b=%d and n=%d and si=%d\n", s, b0, b,
                       n_shelf_cells_to_split[s * m_n_basins + b],
                       add_shelf_instance[s * m_n_basins + b]);
      }
    }
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (cell_type.as_int(i, j) == MASK_FLOATING) {
      int b = basin_mask.as_int(i, j);
      int s = shelf_mask.as_int(i, j);
      if (add_shelf_instance[s * m_n_basins + b] > 0) {
        m_tmp(i, j) = add_shelf_instance[s * m_n_basins + b];
      }
    }
  }

  shelf_mask.copy_from(m_tmp);
}

/*!
 * Compute distance to the grounding line.
 */
void PicoGeometry::compute_distances_gl(const array::Scalar &ocean_mask,
                                        const array::Scalar1 &ice_rises,
                                        bool exclude_ice_rises,
                                        array::Scalar1 &result) {

  array::AccessScope list{ &ice_rises, &ocean_mask, &result };

  result.set(-1);

  // Find the grounding line and the ice front and set result to 1 if ice shelf cell is
  // next to the grounding line, Ice holes within the shelf are treated like ice shelf
  // cells, if exclude_ice_rises is set then ice rises are also treated as ice shelf cells.

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
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

/*!
 * Compute distance to the calving front.
 */
void PicoGeometry::compute_distances_cf(const array::Scalar1 &ocean_mask,
                                        const array::Scalar &ice_rises,
                                        bool exclude_ice_rises,
                                        array::Scalar1 &result) {

  array::AccessScope list{ &ice_rises, &ocean_mask, &result };

  result.set(-1);

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (ice_rises.as_int(i, j) == FLOATING or
          ocean_mask.as_int(i, j) == 1 or
          (exclude_ice_rises and ice_rises.as_int(i, j) == RISE)) {
        // if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

        // label the shelf cells adjacent to the ice front with 1
        auto M = ocean_mask.star(i, j);

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
 *
 * FIXME: replace this with a better algorithm.
 */
void eikonal_equation(array::Scalar1 &mask) {

  assert(mask.stencil_width() > 0);

  auto grid = mask.grid();

  double current_label = 1;
  double continue_loop = 1;
  while (continue_loop != 0) {

    continue_loop = 0;

    for (auto p = grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.as_int(i, j) == 0) {

        auto R = mask.star(i, j);

        if (R.c == 0 and
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

/*!
 * Compute the mask identifying ice shelf "boxes" using distances to the grounding line
 * and the calving front.
 */
void PicoGeometry::compute_box_mask(const array::Scalar &D_gl, const array::Scalar &D_cf,
                                    const array::Scalar &shelf_mask, int max_number_of_boxes,
                                    array::Scalar &result) {

  array::AccessScope list{ &D_gl, &D_cf, &shelf_mask, &result };

  int n_shelves = static_cast<int>(array::max(shelf_mask)) + 1;

  std::vector<double> GL_distance_max(n_shelves, 0.0);
  std::vector<double> GL_distance_max1(n_shelves, 0.0);
  std::vector<double> CF_distance_max(n_shelves, 0.0);
  std::vector<double> CF_distance_max1(n_shelves, 0.0);

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);
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
  GlobalMax(m_grid->com, GL_distance_max.data(), GL_distance_max1.data(), n_shelves);
  GlobalMax(m_grid->com, CF_distance_max.data(), CF_distance_max1.data(), n_shelves);
  // copy data
  GL_distance_max = GL_distance_max1;
  CF_distance_max = CF_distance_max1;

  double GL_distance_ref = *std::max_element(GL_distance_max.begin(), GL_distance_max.end());

  // compute the number of boxes in each shelf

  std::vector<int> n_boxes(n_shelves, 0);
  const int n_min   = 1;
  const double zeta = 0.5;

  for (int k = 0; k < n_shelves; ++k) {
    n_boxes[k] = n_min + round(pow((GL_distance_max[k] / GL_distance_ref), zeta) * (max_number_of_boxes - n_min));

    n_boxes[k] = std::min(n_boxes[k], max_number_of_boxes);
  }

  result.set(0.0);

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    int d_gl = D_gl.as_int(i, j);
    int d_cf = D_cf.as_int(i, j);

    if (shelf_mask.as_int(i, j) > 0 and d_gl > 0.0 and d_cf > 0.0) {
      int shelf_id = shelf_mask.as_int(i, j);
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
