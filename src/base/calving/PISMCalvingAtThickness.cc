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

#include "PISMCalvingAtThickness.hh"
#include "Mask.hh"

PISMCalvingAtThickness::PISMCalvingAtThickness(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent(g, conf) {
  m_calving_threshold = config.get("calving_at_thickness");

  PetscErrorCode ierr = m_old_mask.create(grid, "old_mask", true, 1);
  if (ierr != 0) {
    PetscPrintf(grid.com,
                "PISM ERROR: memory allocation failed (PISMCalvingAtThickness constructor.\n");
    PISMEnd();
  }
}

PISMCalvingAtThickness::~PISMCalvingAtThickness() {
  // empty
}


PetscErrorCode PISMCalvingAtThickness::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 'calving at a threshold thickness' mechanism...\n"
                    "  thickness threshold: %3.3f meters\n", m_calving_threshold); CHKERRQ(ierr);

  return 0;
}

/**
 * Updates ice cover mask and the ice thickness according to the
 * calving rule removing ice at the shelf front that is thinner than a
 * given threshold.
 *
 * @param[in,out] pism_mask ice cover mask
 * @param[in,out] ice_thickness ice thickness
 *
 * @return 0 on success
 */
PetscErrorCode PISMCalvingAtThickness::update(IceModelVec2Int &pism_mask,
                                              IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;
  MaskQuery M(m_old_mask);

  // this call fills ghosts of m_old_mask
  ierr = pism_mask.copy_to(m_old_mask); CHKERRQ(ierr);

  ierr = pism_mask.begin_access();     CHKERRQ(ierr);
  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = m_old_mask.begin_access();    CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.floating_ice(i, j)           &&
          M.next_to_ice_free_ocean(i, j) &&
          ice_thickness(i, j) < m_calving_threshold) {
        ice_thickness(i, j) = 0.0;
        pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
      }
    }
  }

  ierr = m_old_mask.end_access();    CHKERRQ(ierr);
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = pism_mask.end_access();     CHKERRQ(ierr);

  ierr = pism_mask.update_ghosts();     CHKERRQ(ierr);
  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


void PISMCalvingAtThickness::add_vars_to_output(string /*keyword*/,
                                                set<string> &/*result*/) {
  // empty
}

PetscErrorCode PISMCalvingAtThickness::define_variables(set<string> /*vars*/,
                                                        const PIO &/*nc*/,
                                                        PISM_IO_Type /*nctype*/) {
  // empty
  return 0;
}

PetscErrorCode PISMCalvingAtThickness::write_variables(set<string> /*vars*/,
                                                       const PIO& /*nc*/) {
  // empty
  return 0;
}

