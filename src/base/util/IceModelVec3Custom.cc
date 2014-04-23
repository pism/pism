/* Copyright (C) 2013, 2014 PISM Authors
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

#include "IceModelVec3Custom.hh"
#include <assert.h>

namespace pism {

IceModelVec3Custom::IceModelVec3Custom()
{
  // empty
}

IceModelVec3Custom::~IceModelVec3Custom()
{
  // empty
}

/** 
 * Allocate storage and set metadata.
 *
 * @param mygrid grid to use
 * @param short_name name of the NetCDF variable
 * @param z_name name of the NetCDF dimension and variable corresponding to the third dimension
 * @param my_zlevels "vertical" levels (values of z)
 * @param z_attrs attributes of the "z" coordinate variable
 *
 * @return 0 on success
 */

PetscErrorCode IceModelVec3Custom::create(IceGrid &mygrid,
                                          std::string short_name,
                                          std::string z_name,
                                          std::vector<double> my_zlevels,
                                          std::map<std::string, std::string> z_attrs) {
  PetscErrorCode ierr;
  assert(v == NULL);

  m_has_ghosts = false;
  grid         = &mygrid;
  m_name       = short_name;
  zlevels      = my_zlevels;
  m_n_levels   = zlevels.size();

  m_da_stencil_width = 1;

  ierr = grid->get_dm(this->m_n_levels, this->m_da_stencil_width, m_da); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(m_da, &v); CHKERRQ(ierr);

  m_dof = 1;

  m_metadata.resize(m_dof, NCSpatialVariable(grid->get_unit_system()));
  m_metadata[0].init_3d(m_name, mygrid, zlevels);
  m_metadata[0].get_z().set_name(z_name);

  std::map<std::string, std::string>::const_iterator j = z_attrs.begin();
  while (j != z_attrs.end()) {
    if (j->first == "units") {
      m_metadata[0].get_z().set_units(j->second);
    } else if (j->first == "glaciological_units") {
      m_metadata[0].get_z().set_glaciological_units(j->second);
    } else {
      m_metadata[0].get_z().set_string(j->first, j->second);
    }
    ++j;
  }

  return 0;
}

} // end of namespace pism
