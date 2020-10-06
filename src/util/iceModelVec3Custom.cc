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

#include "iceModelVec3Custom.hh"

#include "IceGrid.hh"
#include "error_handling.hh"
#include "pism/util/Context.hh"
#include "pism/util/IceModelVec_impl.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

/**
 * Allocate storage and set metadata.
 *
 * @param grid grid to use
 * @param name name of the NetCDF variable
 * @param z_name name of the NetCDF dimension and variable corresponding to the third dimension
 * @param z_levels "vertical" levels (values of z)
 * @param z_attrs attributes of the "z" coordinate variable
 */
IceModelVec3Custom::IceModelVec3Custom(IceGrid::ConstPtr grid,
                                       const std::string &name,
                                       const std::string &z_name,
                                       const std::vector<double> &z_levels,
                                       const std::map<std::string, std::string> &z_attrs) {
  m_impl->grid = grid;
  m_impl->name = name;
  m_impl->ghosted = false;
  m_impl->zlevels = z_levels;
  m_impl->da_stencil_width = 1;
  m_impl->dof = 1;

  m_impl->da = m_impl->grid->get_dm(m_impl->zlevels.size(), m_impl->da_stencil_width);

  PetscErrorCode ierr = DMCreateGlobalVector(*m_impl->da, m_impl->v.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");

  m_impl->metadata.push_back(SpatialVariableMetadata(m_impl->grid->ctx()->unit_system(),
                                                     m_impl->name, m_impl->zlevels));
  m_impl->metadata[0].get_z().set_name(z_name);

  for (auto z_attr : z_attrs) {
    m_impl->metadata[0].get_z().set_string(z_attr.first, z_attr.second);
  }
}

IceModelVec3Custom::~IceModelVec3Custom() {
  // empty
}

} // end of namespace pism
