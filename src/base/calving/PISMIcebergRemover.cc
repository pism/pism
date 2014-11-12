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

#include "PISMIcebergRemover.hh"
#include "connected_components.hh"
#include "Mask.hh"
#include "PISMVars.hh"
#include "error_handling.hh"

namespace pism {

IcebergRemover::IcebergRemover(IceGrid &g, const Config &conf)
  : Component(g, conf) {

  PetscErrorCode ierr = allocate();
  if (ierr != 0) {
    throw std::runtime_error("IcebergRemover allocation failed");
  }
}

IcebergRemover::~IcebergRemover() {
  PetscErrorCode ierr = deallocate(); CHKERRCONTINUE(ierr);
}

PetscErrorCode IcebergRemover::init(Vars &vars) {
  m_bcflag = vars.get_2d_mask("bcflag");
  return 0;
}

/**
 * Use PISM's ice cover mask to update ice thickness, removing "icebergs".
 *
 * @param[in,out] pism_mask PISM's ice cover mask
 * @param[in,out] ice_thickness ice thickness
 */
PetscErrorCode IcebergRemover::update(IceModelVec2Int &pism_mask,
                                      IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;
  const int
    mask_grounded_ice = 1,
    mask_floating_ice = 2;
  MaskQuery M(pism_mask);

  // prepare the mask that will be handed to the connected component
  // labeling code:
  {
    ierr = m_iceberg_mask.set(0.0); CHKERRQ(ierr);

    IceModelVec::AccessList list;
    list.add(pism_mask);
    list.add(m_iceberg_mask);

    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (M.grounded_ice(i,j) == true) {
        m_iceberg_mask(i,j) = mask_grounded_ice;
      } else if (M.floating_ice(i,j) == true) {
        m_iceberg_mask(i,j) = mask_floating_ice;
      }
    }

    // Mark icy SSA Dirichlet B.C. cells as "grounded" because we
    // don't want them removed.
    if (m_bcflag) {
      list.add(*m_bcflag);

      for (Points p(grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if ((*m_bcflag)(i,j) > 0.5 && M.icy(i,j)) {
          m_iceberg_mask(i,j) = mask_grounded_ice;
        }
      }
    }
  }

  // identify icebergs using serial code on processor 0:
  {
    ierr = m_iceberg_mask.put_on_proc0(m_mask_p0); CHKERRQ(ierr);

    if (grid.rank == 0) {
      double *mask;
      ierr = VecGetArray(m_mask_p0, &mask);
      PISM_PETSC_CHK(ierr, "VecGetArray");

      cc(mask, grid.Mx, grid.My, true, mask_grounded_ice);

      ierr = VecRestoreArray(m_mask_p0, &mask);
      PISM_PETSC_CHK(ierr, "VecRestoreArray");
    }

    ierr = m_iceberg_mask.get_from_proc0(m_mask_p0); CHKERRQ(ierr);
  }

  // correct ice thickness and the cell type mask using the resulting
  // "iceberg" mask:
  {
    IceModelVec::AccessList list;
    list.add(ice_thickness);
    list.add(pism_mask);
    list.add(m_iceberg_mask);

    if (m_bcflag != NULL) {
      // if SSA Dirichlet B.C. are in use, do not modify mask and ice
      // thickness at Dirichlet B.C. locations
      list.add(*m_bcflag);

      for (Points p(grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (m_iceberg_mask(i,j) > 0.5 && (*m_bcflag)(i,j) < 0.5) {
          ice_thickness(i,j) = 0.0;
          pism_mask(i,j)     = MASK_ICE_FREE_OCEAN;
        }
      }
    } else {

      for (Points p(grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (m_iceberg_mask(i,j) > 0.5) {
          ice_thickness(i,j) = 0.0;
          pism_mask(i,j)     = MASK_ICE_FREE_OCEAN;
        }
      }
    }
  }

  // update ghosts of the mask and the ice thickness (then surface
  // elevation can be updated redundantly)
  ierr = pism_mask.update_ghosts(); CHKERRQ(ierr);
  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IcebergRemover::allocate() {
  PetscErrorCode ierr;

  ierr = m_iceberg_mask.create(grid, "iceberg_mask", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_iceberg_mask.allocate_proc0_copy(m_mask_p0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IcebergRemover::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(&m_mask_p0); CHKERRCONTINUE(ierr);

  return 0;
}

void IcebergRemover::add_vars_to_output(const std::string &, std::set<std::string> &) {
  // empty
}

PetscErrorCode IcebergRemover::define_variables(const std::set<std::string> &, const PIO &,
                                                IO_Type) {
  // empty
  return 0;
}

PetscErrorCode IcebergRemover::write_variables(const std::set<std::string> &, const PIO&) {
  // empty
  return 0;
}

} // end of namespace pism
