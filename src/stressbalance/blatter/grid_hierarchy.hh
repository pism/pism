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

/*!
 * A structure wrapping outputs of DMDAGetInfo() and DMDAGetOwnershipRanges().
 */
class DMInfo {
public:
  DMInfo(DM dm);
  DMInfo transpose() const;

  PetscInt Mx, My, Mz;
  PetscInt mx, my, mz;
  PetscInt dims, dof, stencil_width;
  DMDAStencilType stencil_type;
  const PetscInt *lx, *ly, *lz;
  DMBoundaryType bx, by, bz;
private:
  PetscErrorCode setup(DM dm);
};

struct GridInfo {
  double x_min, x_max;
  double y_min, y_max;
  // minimum thickness (used to compute node types)
  double min_thickness;
  // number of 2D parameters on each grid level
  int dof;

  double dx(int mx) const {
    return (x_max - x_min) / (mx - 1);
  }

  double dy(int my) const {
    return (y_max - y_min) / (my - 1);
  }

  double x(double dx, int k) const {
    return x_min + dx * k;
  }

  double y(double dy, int k) const {
    return y_min + dy * k;
  }
};


double grid_z(double b, double H, int Mz, int k);

int grid_padding(int N, int coarsening_factor, int n_levels);

DMDALocalInfo grid_transpose(const DMDALocalInfo &input);

PetscErrorCode setup_level(DM dm, const GridInfo &grid_info);

PetscErrorCode create_restriction(DM fine, DM coarse, const char *dm_name);

PetscErrorCode restrict_data(DM fine, DM coarse, const char *dm_name);

} // end of namespace pism

#endif /* PISM_GRID_HIERARCHY_H */
