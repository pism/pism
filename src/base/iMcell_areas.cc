// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#include "iceModel.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/projection.hh"

namespace pism {

void IceModel::compute_cell_areas() {

  std::string projection = m_grid->get_mapping_info().proj4;

  if (m_config->get_boolean("grid.correct_cell_areas") and
      not projection.empty()) {

    m_log->message(2,
                   "* Computing cell areas using projection parameters (%s)...\n",
                   projection.c_str());

    ::pism::compute_cell_areas(projection, m_cell_area);

    m_log->message(2,
                   "* Computing longitude and latitude using projection parameters (%s)...\n",
                   projection.c_str());

    compute_longitude(projection, m_longitude);
    compute_latitude(projection, m_latitude);
  } else {
    m_log->message(2,
                   "* Computing cell areas using grid spacing (dx = %f m, dy = %f m)...\n",
                   m_grid->dx(), m_grid->dy());

    m_cell_area.set(m_grid->dx() * m_grid->dy());
  }
}

} // end of namespace pism
