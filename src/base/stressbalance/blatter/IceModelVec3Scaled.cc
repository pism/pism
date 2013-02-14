/* Copyright (C) 2013 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
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

#include "IceModelVec3Scaled.hh"

PetscErrorCode IceModelVec3Scaled::create(IceGrid &mygrid, const char my_short_name[], bool local,
					  int myMbz, int stencil_width) {
  PetscErrorCode ierr;
  grid = &mygrid;

  if (!utIsInit()) {
    SETERRQ(grid->com, 1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(grid->com, 2,"IceModelVec3Scaled with name='%s' already allocated\n",name.c_str());
  }

  name = my_short_name;

  n_levels = myMbz;
  zlevels.resize(n_levels);
  double dz = 1.0 / (myMbz - 1);
  for (int i = 0; i < n_levels; ++i)
    zlevels[i] = i * dz;
  zlevels.back() = 1.0;

  da_stencil_width = stencil_width;
  ierr = grid->get_dm(this->n_levels, this->da_stencil_width, da); CHKERRQ(ierr);

  localp = local;
  if (local) {
    ierr = DMCreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  vars[0].init_3d(name, mygrid, zlevels);
  vars[0].dimensions["z"] = "z_scaled";

  map<string,string> &attrs = vars[0].z_attrs;
  attrs["axis"]          = "Z";
  attrs["long_name"]     = "scaled Z-coordinate in the ice (z_base=0, z_surface=1)";
  attrs["units"]         = "1";
  attrs["positive"]      = "up";

  return 0;
}

