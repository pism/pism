/* Copyright (C) 2018, 2020, 2021, 2022, 2023 PISM Authors
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

#ifndef PICOGEOMETRY_H
#define PICOGEOMETRY_H

#include <map>
#include <set>

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace ocean {

void eikonal_equation(array::Scalar1 &mask);

/*!
 * This class isolates geometric computations performed by the PICO ocean model.
 */
class PicoGeometry : public Component {
public:
  PicoGeometry(std::shared_ptr<const Grid> grid);
  virtual ~PicoGeometry() = default;

  void init();
  void update(const array::Scalar &bed_elevation,
              const array::CellType1 &cell_type);

  const array::Scalar &continental_shelf_mask() const;
  const array::Scalar &box_mask() const;
  const array::Scalar &ice_shelf_mask() const;
  const array::Scalar &ice_rise_mask() const;
  const array::Scalar &basin_mask() const;

  enum IceRiseMask { OCEAN = 0, RISE = 1, CONTINENTAL = 2, FLOATING = 3 };

private:
  void compute_ice_rises(const array::CellType &cell_type,
                         bool exclude_ice_rises, array::Scalar &result);
  void compute_lakes(const array::CellType &cell_type, array::Scalar &result);
  void compute_ocean_mask(const array::CellType &cell_type, array::Scalar &result);
  void compute_continental_shelf_mask(const array::Scalar &bed_elevation,
                                      const array::Scalar &ice_rise_mask,
                                      double bed_elevation_threshold,
                                      array::Scalar &result);
  void compute_ice_shelf_mask(const array::Scalar &ice_rise_mask,
                              const array::Scalar &lake_mask,
                              array::Scalar &result);

  std::vector<std::set<int> > basin_neighbors(const array::CellType1 &cell_type,
                                               const array::Scalar1 &basin_mask);

  void identify_calving_front_connection(const array::CellType1 &cell_type,
                                         const array::Scalar &basin_mask,
                                         const array::Scalar &shelf_mask,
                                         int n_shelves,
                                         std::vector<int> &most_shelf_cells_in_basin,
                                         std::vector<int> &cfs_in_basins_per_shelf);

  void split_ice_shelves(const array::CellType &cell_type,
                         const array::Scalar &basin_mask,
                         const std::vector<std::set<int> > &basin_neighbors,
                         const std::vector<int> &most_shelf_cells_in_basin,
                         const std::vector<int> &cfs_in_basins_per_shelf,
                         int n_shelves,
                         array::Scalar &shelf_mask);
 
  void compute_distances_cf(const array::Scalar1 &ocean_mask,
                            const array::Scalar &ice_rises,
                            bool exclude_ice_rises,
                            array::Scalar1 &result);

  void compute_distances_gl(const array::Scalar &ocean_mask,
                            const array::Scalar1 &ice_rises,
                            bool exclude_ice_rises,
                            array::Scalar1 &result);

  void compute_box_mask(const array::Scalar &D_gl,
                        const array::Scalar &D_cf,
                        const array::Scalar &shelf_mask,
                        int max_number_of_boxes,
                        array::Scalar &result);

  void relabel_by_size(array::Scalar &mask);

  // storage for outputs
  array::Scalar m_continental_shelf;
  array::Scalar m_boxes;
  array::Scalar m_ice_shelves;
  array::Scalar1 m_basin_mask;

  // storage for intermediate fields
  array::Scalar1 m_distance_gl;
  array::Scalar1 m_distance_cf;
  array::Scalar1 m_ocean_mask;
  array::Scalar m_lake_mask;
  array::Scalar1 m_ice_rises;

  // temporary storage
  array::Scalar m_tmp;
  std::shared_ptr<petsc::Vec> m_tmp_p0;

  int m_n_basins;
  std::vector<std::set<int> > m_basin_neighbors;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* PICOGEOMETRY_H */
