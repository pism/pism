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

#include "PISMIcebergRemover.hh"
#include "connected_components.hh"
#include "Mask.hh"

PISMIcebergRemover::PISMIcebergRemover(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent(g, conf) {

  PetscErrorCode ierr = allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: failed to allocate PISMIcebergRemover.\n");
    PISMEnd();
  }
}

PISMIcebergRemover::~PISMIcebergRemover() {
  PetscErrorCode ierr = deallocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: failed to deallocate PISMIcebergRemover.\n");
    PISMEnd();
  }
}

PetscErrorCode PISMIcebergRemover::init(PISMVars &/*vars*/) {
  // do nothing
  return 0;
}

/**
 * Use PISM's ice cover mask to update ice thickness, removing "icebergs".
 *
 * @param[in,out] pism_mask PISM's ice cover mask
 * @param[in,out] ice_thickness ice thickness
 */
PetscErrorCode PISMIcebergRemover::update(IceModelVec2Int &pism_mask,
                                          IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;
  PetscScalar **iceberg_mask;
  const int
    mask_grounded_ice = 1,
    mask_floating_ice = 2;
  MaskQuery M(pism_mask);

  // prepare the mask that will be handed to the connected component
  // labeling code:
  ierr = VecSet(m_g2, 0.0); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(m_da2, m_g2, &iceberg_mask); CHKERRQ(ierr);
  ierr = pism_mask.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.grounded_ice(i,j) == true)
        iceberg_mask[i][j] = mask_grounded_ice;
      else if (M.floating_ice(i,j) == true)
        iceberg_mask[i][j] = mask_floating_ice;
    }
  }
  ierr = pism_mask.end_access(); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(m_da2, m_g2, &iceberg_mask); CHKERRQ(ierr);

  ierr = transfer_to_proc0(); CHKERRQ(ierr);

  if (grid.rank == 0) {
    PetscScalar *mask;
    ierr = VecGetArray(m_mask_p0, &mask); CHKERRQ(ierr);

    cc(mask, grid.Mx, grid.My, true, mask_grounded_ice);

    ierr = VecRestoreArray(m_mask_p0, &mask); CHKERRQ(ierr);
  }

  ierr = transfer_from_proc0(); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(m_da2, m_g2, &iceberg_mask); CHKERRQ(ierr);
  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = pism_mask.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (iceberg_mask[i][j] > 0.5) {
        ice_thickness(i,j) = 0.0;
        pism_mask(i,j)     = MASK_ICE_FREE_OCEAN;
      }
    }
  }
  ierr = pism_mask.end_access(); CHKERRQ(ierr);
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(m_da2, m_g2, &iceberg_mask); CHKERRQ(ierr);

  // update ghosts of the mask and the ice thickness (then surface
  // elevation can be updated redundantly)
  ierr = pism_mask.update_ghosts(); CHKERRQ(ierr);
  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMIcebergRemover::allocate() {
  PetscErrorCode ierr;

  ierr = grid.get_dm(1,         // dof
                     1,         // stencil width
                     m_da2); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(m_da2, &m_g2); CHKERRQ(ierr);

  // We want a global Vec but reordered in the natural ordering so
  // when it is scattered to proc zero it is not all messed up
  ierr = DMDACreateNaturalVector(m_da2, &m_g2natural); CHKERRQ(ierr);

  // Get scatter context *and* allocate mask on processor 0:
  ierr = VecScatterCreateToZero(m_g2natural, &m_scatter, &m_mask_p0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMIcebergRemover::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(&m_g2); CHKERRQ(ierr);
  ierr = VecDestroy(&m_g2natural); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&m_scatter); CHKERRQ(ierr);
  ierr = VecDestroy(&m_mask_p0); CHKERRQ(ierr);

  return 0;
}

/**
 * Transfer the m_g2 data member to m_mask_p0 on rank 0.
 */
PetscErrorCode PISMIcebergRemover::transfer_to_proc0() {
  PetscErrorCode ierr;

  ierr = DMDAGlobalToNaturalBegin(m_da2, m_g2, INSERT_VALUES, m_g2natural); CHKERRQ(ierr);
  ierr =   DMDAGlobalToNaturalEnd(m_da2, m_g2, INSERT_VALUES, m_g2natural); CHKERRQ(ierr);

  ierr = VecScatterBegin(m_scatter, m_g2natural, m_mask_p0,
                         INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr =   VecScatterEnd(m_scatter, m_g2natural, m_mask_p0,
                         INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  return 0;
}


/**
 * Transfer the m_mask_p0 data member from rank 0 to m_g2 (distributed).
 */
PetscErrorCode PISMIcebergRemover::transfer_from_proc0() {
  PetscErrorCode ierr;

  ierr = VecScatterBegin(m_scatter, m_mask_p0, m_g2natural,
                         INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
  ierr =   VecScatterEnd(m_scatter, m_mask_p0, m_g2natural,
                         INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = DMDANaturalToGlobalBegin(m_da2, m_g2natural, INSERT_VALUES, m_g2); CHKERRQ(ierr);
  ierr =   DMDANaturalToGlobalEnd(m_da2, m_g2natural, INSERT_VALUES, m_g2); CHKERRQ(ierr);

  return 0;
}

void PISMIcebergRemover::add_vars_to_output(string, set<string> &) {
  // empty
}

PetscErrorCode PISMIcebergRemover::define_variables(set<string>, const PIO &,
                                                    PISM_IO_Type) {
  // empty
  return 0;
}

PetscErrorCode PISMIcebergRemover::write_variables(set<string>, const PIO& ) {
  // empty
  return 0;
}
