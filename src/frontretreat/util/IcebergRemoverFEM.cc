/* Copyright (C) 2021, 2022, 2023 PISM Authors
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

#include <cassert>

#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/connected_components.hh"

#include "pism/frontretreat/util/IcebergRemoverFEM.hh"

#include "pism/util/fem/Element.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/Mask.hh"
#include "pism/util/interpolation.hh"

namespace pism {
namespace calving {

IcebergRemoverFEM::IcebergRemoverFEM(std::shared_ptr<const Grid> grid)
  : IcebergRemover(grid),
    m_mask(grid, "temporary_mask") {
  m_mask.set_interpolation_type(NEAREST);
}

/*! Remove "icebergs" using the finite element notion of connectivity: two elements are
 *  connected if they share a boundary.
 *
 * 1. Loop over elements and create a mask that will be used to determine connectivity
 *    between elements.
 *
 * - an element is a "grounded ice element" if all nodes are icy and are either grounded
 *   or belong to the set of Dirichlet nodes
 *
 *    - an element is "floating ice" if all nodes are icy and at least one node is
 *      "floating ice"
 *
 *    - all other elements are ice-free
 *
 * 2. Label connected components, identifying "icebergs".
 *
 * Once "iceberg" elements are labeled we need to remove *nodes* that belong to icebergs
 * but *do not* belong to any elements connected to grounded ice.
 *
 * 3. Create a mask filled with zeros. Loop over elements and add 1 to nodes of all
 *    "iceberg" elements. Add -1 to all nodes of "grounded" elements.
 *
 * 4. Now loop over all nodes and remove nodes with positive mask values.
 *
 */
void IcebergRemoverFEM::update_impl(const array::Scalar &bc_mask,
                                    array::CellType1 &cell_type,
                                    array::Scalar &ice_thickness) {
  const int
    mask_grounded_ice = 1,
    mask_floating_ice = 2;

  int bc_mask_nodal[fem::q1::n_chi];
  int cell_type_nodal[fem::q1::n_chi];

  assert(bc_mask.stencil_width() >= 1);
  assert(cell_type.stencil_width() >= 1);

  array::AccessScope list{&bc_mask, &cell_type, &m_iceberg_mask};

  fem::Q1Element2 element(*m_grid, fem::Q1Quadrature1());

  // prepare the iceberg mask: the value at (i, j) describes an *element* with (i,j) as a
  // lower left corner
  {
    // loop over all nodes in a local sub-domain
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      element.reset(i, j);
      // the following two calls use ghost values
      element.nodal_values(bc_mask, bc_mask_nodal);
      element.nodal_values(cell_type, cell_type_nodal);

      // check if all nodes are icy
      bool icy = true;
      for (int n = 0; icy and n < fem::q1::n_chi; ++n) {
        icy &= mask::icy(cell_type_nodal[n]);
      }

      if (icy) {
        // This is an icy element: check if all nodes are grounded or are a part of the
        // set of Dirichlet nodes
        bool grounded = true;
        for (int n = 0; grounded and n < fem::q1::n_chi; ++n) {
          grounded &= (mask::grounded(cell_type_nodal[n]) or bc_mask_nodal[n] == 1);
        }

        m_iceberg_mask(i, j) = grounded ? mask_grounded_ice : mask_floating_ice;
      } else {
        // This is an ice-free element.
        m_iceberg_mask(i, j) = 0;
      }
    } // end of the loop over local nodes
  } // end of the block preparing the mask

  // Identify icebergs using serial code on processor 0:
  {
    m_iceberg_mask.put_on_proc0(*m_mask_p0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        petsc::VecArray mask_p0(*m_mask_p0);
        label_connected_components(mask_p0.get(), m_grid->My(), m_grid->Mx(),
                                   true, mask_grounded_ice);
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();

    m_iceberg_mask.get_from_proc0(*m_mask_p0);
    // note: this will update ghosts of m_iceberg_mask
  }

  // create a mask indicating if a *node* should be removed
  {
    DMDALocalInfo info;
    {
      auto da = m_grid->get_dm(1, 0);  // dof = 1, stencil_width = 0
      PetscErrorCode ierr = DMDAGetLocalInfo(*da, &info);
      if (ierr != 0) {
        throw std::runtime_error("Failed to get DMDA info");
      }
    }

    m_mask.set(0);
    list.add(m_mask);
    double **M = m_mask.array();

    double mask_iceberg[] = {1.0, 1.0, 1.0, 1.0};
    double mask_grounded[] = {-1.0, -1.0, -1.0, -1.0};

    // loop over all the elements that have at least one owned node
    for (int j = info.gys; j < info.gys + info.gym - 1; j++) {
      for (int i = info.gxs; i < info.gxs + info.gxm - 1; i++) {
        element.reset(i, j);

        // the following two calls use ghost values
        element.nodal_values(bc_mask, bc_mask_nodal);
        element.nodal_values(cell_type, cell_type_nodal);

        // check if all nodes are icy
        bool icy = true;
        for (int n = 0; icy and n < fem::q1::n_chi; ++n) {
          icy &= mask::icy(cell_type_nodal[n]);
        }

        if (icy) {
          // check if all nodes are grounded or are a part of the set of Dirichlet nodes
          bool grounded = true;
          for (int n = 0; grounded and n < fem::q1::n_chi; ++n) {
            grounded &= (mask::grounded(cell_type_nodal[n]) or bc_mask_nodal[n] == 1);
          }

          if (m_iceberg_mask(i, j) == 1) {
            // this is an iceberg element
            element.add_contribution(mask_iceberg, M);
          } else {
            element.add_contribution(mask_grounded, M);
          }
        }
      }
    } // end of the loop over elements
  } // end of the block identifying nodes to remove

  // loop over all *nodes* and modify ice thickness and mask
  {
    list.add(ice_thickness);

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_mask(i, j) > 0) {
        ice_thickness(i,j) = 0.0;
        cell_type(i,j)          = MASK_ICE_FREE_OCEAN;
      }
    }
  }

  // update ghosts of the mask and the ice thickness (then surface
  // elevation can be updated redundantly)
  cell_type.update_ghosts();
  ice_thickness.update_ghosts();
}

} // end of namespace calving
} // end of namespace pism
