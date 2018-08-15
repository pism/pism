// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "ExpandSL.hh"

namespace pism {
namespace ocean {
namespace sea_level {

ExpandSL::ExpandSL(IceGrid::ConstPtr grid, std::shared_ptr<SeaLevel> in)
  : SeaLevel(grid, in) {
  // empty
}

ExpandSL::~ExpandSL() {
  // empty
}

void ExpandSL::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);
}

void ExpandSL::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);
  m_sea_level.copy_from(m_input_model->elevation());

  IceModelVec2S sea_level_old(m_grid, "sea_level_backup", WITH_GHOSTS, 1);
  sea_level_old.copy_from(m_sea_level);
  sea_level_old.update_ghosts();

  const Direction dirs[] = { North, East, South, West };

  IceModelVec::AccessList list{ &m_sea_level, &sea_level_old };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ( sea_level_old(i, j) != m_fill_value ) {
      StarStencil<double> sl_star = sea_level_old.star(i, j);
      for (int n = 0; n < 4; ++n) {
        const Direction direction = dirs[n];
        const double Level = sl_star[direction];
        if ( Level == m_fill_value ) {
          m_sea_level(i, j) = m_fill_value;
          break;
        }
      } // end loop neighbors
    } // end if
  } // end loop grid

  m_sea_level.update_ghosts();
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
