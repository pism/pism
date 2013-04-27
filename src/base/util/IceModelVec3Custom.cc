/* Copyright (C) 2013 PISM Authors
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
                                          string short_name,
                                          string z_name,
                                          vector<double> my_zlevels,
                                          map<string, string> z_attrs) {
  PetscErrorCode ierr;

  if (v != PETSC_NULL) {
    SETERRQ1(grid->com,
             2,
             "IceModelVec3Custom with name='%s' already allocated\n",
             name.c_str());
  }

  localp   = false;
  grid     = &mygrid;
  name     = short_name;
  zlevels  = my_zlevels;
  n_levels = zlevels.size();

  da_stencil_width = 1;

  ierr = grid->get_dm(this->n_levels, this->da_stencil_width, da); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da, &v); CHKERRQ(ierr);

  vars[0].init_3d(name, mygrid, zlevels);
  vars[0].dimensions["z"] = z_name;
  vars[0].z_attrs         = z_attrs;

  return 0;
}
