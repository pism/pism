/* Copyright (C) 2018, 2020, 2021, 2022 PISM Authors
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
#include "pism/util/IceModelVec2S.hh"

namespace pism {

class IceModelVec2S;

namespace ocean {

void eikonal_equation(IceModelVec2S &mask);

/*!
 * This class isolates geometric computations performed by the PICO ocean model.
 */
class PicoGeometry : public Component {
public:
  PicoGeometry(IceGrid::ConstPtr grid);
  virtual ~PicoGeometry() = default;

  void init();
  void update(const IceModelVec2S &bed_elevation,
              const CellTypeArray1 &cell_type);

  const IceModelVec2S &continental_shelf_mask() const;
  const IceModelVec2S &box_mask() const;
  const IceModelVec2S &ice_shelf_mask() const;
  const IceModelVec2S &ice_rise_mask() const;
  const IceModelVec2S &basin_mask() const;

  enum IceRiseMask { OCEAN = 0, RISE = 1, CONTINENTAL = 2, FLOATING = 3 };

private:
  void compute_ice_rises(const CellTypeArray0 &cell_type,
                         bool exclude_ice_rises, IceModelVec2S &result);
  void compute_lakes(const CellTypeArray0 &cell_type, IceModelVec2S &result);
  void compute_ocean_mask(const CellTypeArray0 &cell_type, IceModelVec2S &result);
  void compute_continental_shelf_mask(const IceModelVec2S &bed_elevation,
                                      const IceModelVec2S &ice_rise_mask,
                                      double bed_elevation_threshold,
                                      IceModelVec2S &result);
  void compute_ice_shelf_mask(const IceModelVec2S &ice_rise_mask,
                              const IceModelVec2S &lake_mask,
                              IceModelVec2S &result);

  std::map<int,std::set<int> > basin_neighbors(const CellTypeArray1 &cell_type,
                                               const IceModelVec2S &basin_mask);

  void identify_calving_front_connection(const CellTypeArray1 &cell_type,
                                         const IceModelVec2S &basin_mask,
                                         const IceModelVec2S &shelf_mask,
                                         int n_shelves,
                                         std::vector<int> &most_shelf_cells_in_basin,
                                         std::vector<int> &cfs_in_basins_per_shelf);

  void split_ice_shelves(const CellTypeArray0 &cell_type,
                         const IceModelVec2S &basin_mask,
                         const std::map<int, std::set<int> > &basin_neighbors,
                         const std::vector<int> &most_shelf_cells_in_basin,
                         const std::vector<int> &cfs_in_basins_per_shelf,
                         int n_shelves,
                         IceModelVec2S &shelf_mask);
 
  void compute_distances_cf(const IceModelVec2S &ocean_mask,
                            const IceModelVec2S &ice_rises,
                            bool exclude_ice_rises,
                            IceModelVec2S &result);

  void compute_distances_gl(const IceModelVec2S &ocean_mask,
                            const IceModelVec2S &ice_rises,
                            bool exclude_ice_rises,
                            IceModelVec2S &result);

  void compute_box_mask(const IceModelVec2S &D_gl,
                        const IceModelVec2S &D_cf,
                        const IceModelVec2S &shelf_mask,
                        int max_number_of_boxes,
                        IceModelVec2S &result);

  void label_tmp();
  void relabel_by_size(IceModelVec2S &mask);

  // storage for outputs
  IceModelVec2S m_continental_shelf;
  IceModelVec2S m_boxes;
  IceModelVec2S m_ice_shelves;
  Array2SGhosted<1> m_basin_mask;

  // storage for intermediate fields
  Array2SGhosted<1> m_distance_gl;
  Array2SGhosted<1> m_distance_cf;
  Array2SGhosted<1> m_ocean_mask;
  IceModelVec2S m_lake_mask;
  Array2SGhosted<1> m_ice_rises;

  // temporary storage
  IceModelVec2S m_tmp;
  std::shared_ptr<petsc::Vec> m_tmp_p0;

  int m_n_basins;
  std::map<int, std::set<int> > m_basin_neighbors;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* PICOGEOMETRY_H */
