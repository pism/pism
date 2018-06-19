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

#include "pism/util/Vars.hh"
#include "pism/geometry/Geometry.hh"

#include "ExpandLakeMargins.hh"

namespace pism {
namespace ocean {
namespace lake_level {

ExpandLakeMargins::ExpandLakeMargins(IceGrid::ConstPtr grid,
                                     std::shared_ptr<LakeLevel> in)
  : LakeLevel(grid, in) {
  //empty
}


void ExpandLakeMargins::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);
}

void ExpandLakeMargins::update_impl(const Geometry &geometry, double t, double dt) {

  m_input_model->update(geometry, t, dt);
  m_lake_level.copy_from(m_input_model->elevation());

  IceModelVec2S lake_level_old(m_grid, "lake_level_backup", WITH_GHOSTS, 1);
  lake_level_old.copy_from(m_lake_level);
  lake_level_old.update_ghosts();

  const Direction dirs[] = { North, East, South, West };

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list{ &m_lake_level, &lake_level_old };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (not gc.islake(lake_level_old(i, j))) {
      StarStencil<double> ll_star = lake_level_old.star(i, j);
      for (int n = 0; n < 4; ++n) {
        const Direction direction = dirs[n];
        const double Level = ll_star[direction];
        if (gc.islake(Level)) {
          m_lake_level(i, j) = Level;
          break;
        }
      } // end loop neighbors
    } // end "if not lake"
  } // end loop grid

  m_lake_level.update_ghosts();
}


} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
