/* Copyright (C) 2020 PISM Authors
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
#include <petscmat.h>
#include <petscvec.h>

namespace pism {

struct GridInfo {
  // half-width of the domain in the X direction
  double Lx;
  // half-width of the domain in the Y direction
  double Ly;
  // minimum thickness (used to compute node types)
  double min_thickness;
  // number of degrees of freedom in the 2D input Vec
  int dof;
};

PetscErrorCode setup_level(DM dm, const GridInfo &grid_info);

PetscErrorCode create_restriction(DM fine, DM coarse,
                                  const char *dm_name, const char *mat_name);

PetscErrorCode restrict_data(DM fine, DM coarse,
                             const char *dm_name,
                             const char *mat_name,
                             const char *vec_name);

PetscErrorCode restriction_hook(DM fine,
                                Mat mrestrict, Vec rscale, Mat inject,
                                DM coarse, void *ctx);

PetscErrorCode coarsening_hook(DM dm_fine, DM dm_coarse, void *ctx);

} // end of namespace pism

#endif /* PISM_GRID_HIERARCHY_H */
