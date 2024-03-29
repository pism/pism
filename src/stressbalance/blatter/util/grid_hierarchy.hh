/* Copyright (C) 2020, 2021 PISM Authors
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

#ifndef PISM_GRID_HIERARCHY_H
#define PISM_GRID_HIERARCHY_H

#include <petscdmda.h>

namespace pism {

double grid_z(double b, double H, int Mz, int k);

DMDALocalInfo grid_transpose(const DMDALocalInfo &input);

PetscErrorCode setup_level(DM dm, int mg_levels);

PetscErrorCode create_restriction(DM fine, DM coarse, const char *dm_name);

PetscErrorCode restrict_data(DM fine, DM coarse, const char *dm_name);

} // end of namespace pism

#endif /* PISM_GRID_HIERARCHY_H */
