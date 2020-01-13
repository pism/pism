/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2019, 2020 PISM Authors
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

#include <cassert>

#include "iceModelVec3Custom.hh"
#include "ConfigInterface.hh"
#include "IceGrid.hh"
#include "error_handling.hh"

namespace pism {

/**
 * Allocate storage and set metadata.
 *
 * @param grid grid to use
 * @param name name of the NetCDF variable
 * @param z_name name of the NetCDF dimension and variable corresponding to the third dimension
 * @param z_levels "vertical" levels (values of z)
 * @param z_attrs attributes of the "z" coordinate variable
 *
 * @return 0 on success
 */
IceModelVec3Custom::IceModelVec3Custom(IceGrid::ConstPtr grid,
                                       const std::string &name,
                                       const std::string &z_name,
                                       const std::vector<double> &z_levels,
                                       const std::map<std::string, std::string> &z_attrs) {
  m_grid = grid;
  m_name = name;
  m_has_ghosts = false;
  m_zlevels = z_levels;
  m_da_stencil_width = 1;
  m_dof = 1;

  m_da = m_grid->get_dm(this->m_zlevels.size(), this->m_da_stencil_width);

  PetscErrorCode ierr = DMCreateGlobalVector(*m_da, m_v.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");

  m_metadata.push_back(SpatialVariableMetadata(m_grid->ctx()->unit_system(),
                                               m_name, m_zlevels));
  m_metadata[0].get_z().set_name(z_name);

  for (auto z_attr : z_attrs) {
    m_metadata[0].get_z().set_string(z_attr.first, z_attr.second);
  }
}

IceModelVec3Custom::~IceModelVec3Custom() {
  // empty
}

} // end of namespace pism
