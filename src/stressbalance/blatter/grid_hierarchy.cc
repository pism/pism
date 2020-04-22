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

namespace pism {

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
                       &Mx, &My, &Mz, // grid size
                       &mx, &my, &mz, // number of processors in each direction
                       NULL,          // number of degrees of freedom
                       &stencil_width,
                       NULL, NULL, NULL, // types of ghost nodes at the boundary
                       &stencil_type); CHKERRQ(ierr);

    ierr = DMDAGetOwnershipRanges(dm, &lx, &ly, &lz); CHKERRQ(ierr);
  }

  // Create a 2D DMDA and a global Vec, then stash them in dm.
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
    int n_3d_dof = 1;

    ierr = DMDACreate3d(comm,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        stencil_type,
                        Mx, My, Mz,
                        mx, my, mz,
                        n_3d_dof,
                        stencil_width,
                        lx, ly, lz,
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

/*!
 * Restrict 2D and 3D model parameters from a fine grid to a coarse grid.
 *
 * This hook is called every time SNES needs to update coarse-grid data.
 */
PetscErrorCode restriction_hook(DM fine,
                                Mat mrestrict, Vec rscale, Mat inject,
                                DM coarse, void *ctx) {
  // Get rid of warnings about unused arguments
  (void) mrestrict;
  (void) rscale;
  (void) inject;
  (void) ctx;

  PetscErrorCode ierr;
  ierr = restrict_data(fine, coarse, "2D_DM"); CHKERRQ(ierr);
  ierr = restrict_data(fine, coarse, "3D_DM"); CHKERRQ(ierr);

  return 0;
}

/*! \brief Grid coarsening hook.
 *
 * This hook is called *once* when SNES sets up the next coarse level.
 *
 * This hook does three things:
 * - Set up the DM for the newly created coarse level.
 * - Set up the matrix type on the coarsest level to allow using
 *   direct solvers for the coarse problem.
 * - Set up the interpolation matrix that will be used by the
 *   restriction hook to set model parameters on the new coarse level.
 *
 * See restriction_hook().
 */
PetscErrorCode coarsening_hook(DM dm_fine, DM dm_coarse, void *ctx) {
  PetscErrorCode ierr;
  GridInfo *grid_info = (GridInfo*)ctx;
  PetscInt rlevel, clevel;

  ierr = setup_level(dm_coarse, *grid_info); CHKERRQ(ierr);

  ierr = DMGetRefineLevel(dm_coarse, &rlevel); CHKERRQ(ierr);
  ierr = DMGetCoarsenLevel(dm_coarse, &clevel); CHKERRQ(ierr);
  if (rlevel - clevel == 0) {
    ierr = DMSetMatType(dm_coarse, MATAIJ); CHKERRQ(ierr);
  }

  ierr = DMCoarsenHookAdd(dm_coarse, coarsening_hook, restriction_hook, ctx); CHKERRQ(ierr);

  // 2D
  ierr = create_restriction(dm_fine, dm_coarse, "2D_DM"); CHKERRQ(ierr);

  // 3D
  ierr = create_restriction(dm_fine, dm_coarse, "3D_DM"); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
