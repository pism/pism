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

#include "grid_hierarchy.hh"
#include "DataAccess.hh"
#include "pism/util/fem/FEM.hh"
#include "pism/util/node_types.hh"

namespace pism {

/*!
 * x and y coordinates of the nodes
 *
 * @param[in] min minimum coordinate value
 * @param[in] delta grid spacing
 * @param[in] k node index
 */
double grid_xy(double min, double delta, int k) {
  return min + k * delta;
}

/*!
 * z coordinates of the nodes
 *
 * @param[in] b surface elevation of the bottom of the domain
 * @param[in] H domain thickness
 * @param[in] Mz number of grid points in each vertical column
 * @param[in] k node index in the z direction
 */
double grid_z(double b, double H, int Mz, int k) {
  return b + H * k / (Mz - 1.0);
}

/*!
 * Compute the padding needed to allow for `n_levels` of coarsening.
 *
 * @param[in] N number of grid points (nodes)
 * @param[in] factor coarsening factor
 * @param[in] n_levels number of coarsening levels
 *
 * @return padding amount
 */
int grid_padding(int N, int factor, int n_levels) {
  // number of spaces
  int k = N - 1;
  int C = 1;
  for (int n = 0; n < n_levels; ++n) {
    C *= factor;
    k = (k % factor ? k + (factor - (k % factor)): k) / factor;
  }
  return (C * k + 1) - N;
}

/* Transpose a DMDALocalInfo structure to map from PETSc's ordering to PISM's order needed
   to ensure that vertical columns are stored contiguously in RAM.

   (What a pain.)

   Map from PETSc to PISM order:

   | PETSc | PISM |
   |-------+------|
   | x     | z    |
   | y     | x    |
   | z     | y    |

   Assuming that i, j, k indexes correspond to x, y, and z PETSc's indexing order is
   [k][j][i]. After this transpose we have to use [j][i][k].

   Note that this indexing order is compatible with the PETSc-standard indexing for 2D
   Vecs: [j][i].

   All the lines changed to implement this transpose are marked with STORAGE_ORDER: that
   way you can use grep to find them.
 */
DMDALocalInfo grid_transpose(const DMDALocalInfo &input) {
  DMDALocalInfo result = input;

  result.mx = input.my;
  result.my = input.mz;
  result.mz = input.mx;

  result.xs = input.ys;
  result.ys = input.zs;
  result.zs = input.xs;

  result.xm = input.ym;
  result.ym = input.zm;
  result.zm = input.xm;

  result.gxs = input.gys;
  result.gys = input.gzs;
  result.gzs = input.gxs;

  result.gxm = input.gym;
  result.gym = input.gzm;
  result.gzm = input.gxm;

  result.bx = input.by;
  result.by = input.bz;
  result.bz = input.bx;

  return result;
}

/*!
 * Set up storage for 2D and 3D data inputs (DMDAs and Vecs)
 */
PetscErrorCode setup_level(DM dm, const GridInfo &grid_info) {
  PetscErrorCode ierr;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm, &comm); CHKERRQ(ierr);

  // Get grid information
  PetscInt Mx, My, Mz, mx, my, mz, stencil_width;
  DMDAStencilType stencil_type;
  const PetscInt *lx, *ly, *lz;
  {
    ierr = DMDAGetInfo(dm,
                       NULL,          // dimensions
                       &Mz, &Mx, &My, // grid size (STORAGE_ORDER)
                       &mz, &mx, &my, // number of processors in each direction (STORAGE_ORDER)
                       NULL,          // number of degrees of freedom
                       &stencil_width,
                       NULL, NULL, NULL, // types of ghost nodes at the boundary
                       &stencil_type); CHKERRQ(ierr);

    ierr = DMDAGetOwnershipRanges(dm, &lz, &lx, &ly); CHKERRQ(ierr); // STORAGE_ORDER
  }

  // Create a 2D DMDA and a global Vec, then stash them in dm.
  //
  // FIXME: we don't need the 2D DMDA now that we don't coarsen in horizontal directions.
  {
    DM  da;
    Vec parameters;

    ierr = DMDACreate2d(comm,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        stencil_type,
                        Mx, My,
                        mx, my,
                        grid_info.dof,
                        stencil_width,
                        lx, ly,
                        &da); CHKERRQ(ierr);

    ierr = DMSetUp(da); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da, &parameters); CHKERRQ(ierr);

    ierr = PetscObjectCompose((PetscObject)dm, "2D_DM", (PetscObject)da); CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)dm, "2D_DM_data", (PetscObject)parameters); CHKERRQ(ierr);

    ierr = DMDestroy(&da); CHKERRQ(ierr);

    ierr = VecDestroy(&parameters); CHKERRQ(ierr);
  }

  // Create a 3D DMDA and a global Vec, then stash them in dm.
  {
    DM  da;
    Vec parameters;
    int dof = 1;

    ierr = DMDACreate3d(comm,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, // STORAGE_ORDER
                        stencil_type,
                        Mz, Mx, My, // STORAGE_ORDER
                        mz, mx, my, // STORAGE_ORDER
                        dof,
                        stencil_width,
                        lz, lx, ly, // STORAGE_ORDER
                        &da); CHKERRQ(ierr);

    ierr = DMSetUp(da); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da, &parameters); CHKERRQ(ierr);

    ierr = PetscObjectCompose((PetscObject)dm, "3D_DM", (PetscObject)da); CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)dm, "3D_DM_data", (PetscObject)parameters); CHKERRQ(ierr);

    ierr = DMDestroy(&da); CHKERRQ(ierr);

    ierr = VecDestroy(&parameters); CHKERRQ(ierr);
  }

  // get refinement level
  PetscInt level = 0;
  ierr = DMGetCoarsenLevel(dm, &level); CHKERRQ(ierr);

  // report
  {
    double
      Wx = grid_info.x_max - grid_info.x_min,
      Wy = grid_info.y_max - grid_info.y_min;
    ierr = PetscPrintf(comm,
                       "Level %D domain size (m) %8.2g x %8.2g,"
                       " num elements %3d x %3d x %3d (%8d), size (m) %g x %g\n",
                       level, Wx, Wy,
                       Mx, My, Mz, Mx*My*Mz, Wx / (Mx - 1), Wy / (My - 1)); CHKERRQ(ierr);
  }
  return 0;
}

/*! @brief Create the restriction matrix.
 *
 * The result of this call is attached to `dm_fine` under `mat_name`.
 *
 * @param[in] fine DM corresponding to the fine grid
 * @param[in] coarse DM corresponding to the coarse grid
 * @param[in] dm_name name of the DM for 2D or 3D parameters
 * @param[in] mat_name name to use when attaching the restriction matrix to `fine`
 */
PetscErrorCode create_restriction(DM fine, DM coarse, const char *dm_name) {
  PetscErrorCode ierr;
  DM da_fine, da_coarse;
  Mat mat;
  Vec scale;

  std::string
    prefix = dm_name,
    mat_name = prefix + "_restriction",
    vec_name = prefix + "_scaling";

  /* Get the DM for parameters from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, dm_name, (PetscObject *)&da_fine); CHKERRQ(ierr);
  if (!da_fine) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No %s composed with given DMDA", dm_name);
  }

  /* 2. get the DM for parameters from the coarse grid DM */
  ierr = PetscObjectQuery((PetscObject)coarse, dm_name, (PetscObject *)&da_coarse); CHKERRQ(ierr);
  if (!da_coarse) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No %s composed with given DMDA", dm_name);
  }

  /* call DMCreateInterpolation */
  ierr = DMCreateInterpolation(da_coarse, da_fine, &mat, &scale); CHKERRQ(ierr);

  /* attach to the fine grid DM */
  ierr = PetscObjectCompose((PetscObject)fine, vec_name.c_str(), (PetscObject)scale); CHKERRQ(ierr);
  ierr = VecDestroy(&scale); CHKERRQ(ierr);

  ierr = PetscObjectCompose((PetscObject)fine, mat_name.c_str(), (PetscObject)mat); CHKERRQ(ierr);
  ierr = MatDestroy(&mat); CHKERRQ(ierr);

  return 0;
}


/*! @brief Restrict model parameters from the `fine` grid onto the `coarse` grid.
 *
 * This function uses the restriction matrix created by coarsening_hook().
 */
PetscErrorCode restrict_data(DM fine, DM coarse, const char *dm_name) {
  PetscErrorCode ierr;
  Vec X_fine, X_coarse, scaling;
  DM da_fine, da_coarse;
  Mat mat;

  std::string
    prefix = dm_name,
    mat_name = prefix + "_restriction",
    scaling_name = prefix + "_scaling",
    vec_name = prefix + "_data";

  /* get the restriction matrix from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, mat_name.c_str(), (PetscObject *)&mat); CHKERRQ(ierr);
  if (!mat) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Failed to get the restriction matrix");
  }

  /* get the scaling vector from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, scaling_name.c_str(), (PetscObject *)&scaling); CHKERRQ(ierr);
  if (!scaling) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Failed to get the scaling vector");
  }

  /* get the DMDA from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, dm_name, (PetscObject *)&da_fine); CHKERRQ(ierr);
  if (!da_fine) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Failed to get the fine grid DM");
  }

  /* get the storage vector from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, vec_name.c_str(), (PetscObject *)&X_fine); CHKERRQ(ierr);
  if (!X_fine) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Failed to get the fine grid Vec");
  }

  /* get the DMDA from the coarse grid DM */
  ierr = PetscObjectQuery((PetscObject)coarse, dm_name, (PetscObject *)&da_coarse); CHKERRQ(ierr);
  if (!da_coarse) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Failed to get the coarse grid DM");
  }

  /* get the storage vector from the coarse grid DM */
  ierr = PetscObjectQuery((PetscObject)coarse, vec_name.c_str(), (PetscObject *)&X_coarse); CHKERRQ(ierr);
  if (!X_coarse) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Failed to get the coarse grid Vec");
  }

  ierr = MatRestrict(mat, X_fine, X_coarse); CHKERRQ(ierr);

  ierr = VecPointwiseMult(X_coarse, X_coarse, scaling); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
