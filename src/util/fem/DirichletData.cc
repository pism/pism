/* Copyright (C) 2020, 2022, 2023 PISM Authors
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

#include "DirichletData.hh"

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/Context.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace fem {

DirichletData::DirichletData()
  : m_indices(NULL), m_weight(1.0) {
  for (unsigned int k = 0; k < q1::n_chi; ++k) {
    m_indices_e[k] = 0;
  }
}

DirichletData::~DirichletData() {
  finish(NULL);
}

void DirichletData::init(const array::Scalar *indices,
                         const array::Array *values,
                         double weight) {
  m_weight = weight;

  if (indices != NULL) {
    indices->begin_access();
    m_indices = indices;
  }

  if (values != NULL) {
    values->begin_access();
  }
}

void DirichletData::finish(const array::Array *values) {
  if (m_indices != NULL) {
    MPI_Comm com = m_indices->grid()->ctx()->com();
    try {
      m_indices->end_access();
      m_indices = NULL;
    } catch (...) {
      handle_fatal_errors(com);
    }
  }

  if (values != NULL) {
    MPI_Comm com = values->grid()->ctx()->com();
    try {
      values->end_access();
    } catch (...) {
      handle_fatal_errors(com);
    }
  }
}

//! @brief Constrain `element`, i.e. ensure that quadratures do not contribute to Dirichlet nodes by marking corresponding rows and columns as "invalid".
void DirichletData::constrain(Element2 &element) {
  element.nodal_values(m_indices->array(), m_indices_e);
  auto n_chi = element.n_chi();
  for (int k = 0; k < n_chi; k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      // Mark any kind of Dirichlet node as not to be touched
      element.mark_row_invalid(k);
      element.mark_col_invalid(k);
    }
  }
}

// Scalar version

DirichletData_Scalar::DirichletData_Scalar(const array::Scalar *indices,
                                           const array::Scalar *values,
                                           double weight)
  : m_values(values) {
  init(indices, m_values, weight);
}

void DirichletData_Scalar::enforce(const Element2 &element, double* x_nodal) {
  assert(m_values != NULL);

  element.nodal_values(m_indices->array(), m_indices_e);
  for (int k = 0; k < element.n_chi(); k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i = 0, j = 0;
      element.local_to_global(k, i, j);
      x_nodal[k] = (*m_values)(i, j);
    }
  }
}

void DirichletData_Scalar::enforce_homogeneous(const Element2 &element, double* x_nodal) {
  element.nodal_values(m_indices->array(), m_indices_e);
  for (int k = 0; k < element.n_chi(); k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_nodal[k] = 0.0;
    }
  }
}

void DirichletData_Scalar::fix_residual(double const *const *const x_global, double **r_global) {
  assert(m_values != NULL);

  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (auto p = grid.points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i] = m_weight * (x_global[j][i] - (*m_values)(i,j));
    }
  }
}

void DirichletData_Scalar::fix_residual_homogeneous(double **r_global) {
  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (auto p = grid.points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i] = 0.0;
    }
  }
}

void DirichletData_Scalar::fix_jacobian(Mat J) {
  const IceGrid &grid = *m_indices->grid();

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.

  const double identity = m_weight;
  ParallelSection loop(grid.com);
  try {
    for (auto p = grid.points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*m_indices)(i, j) > 0.5) {
        MatStencil row;
        row.j = j;
        row.i = i;
        PetscErrorCode ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, &identity,
                                                         ADD_VALUES);
        PISM_CHK(ierr, "MatSetValuesBlockedStencil"); // this may throw
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

DirichletData_Scalar::~DirichletData_Scalar() {
  finish(m_values);
  m_values = NULL;
}

// Vector version

DirichletData_Vector::DirichletData_Vector(const array::Scalar *indices,
                                           const array::Vector *values,
                                           double weight)
  : m_values(values) {
  init(indices, m_values, weight);
}

void DirichletData_Vector::enforce(const Element2 &element, Vector2d* x_nodal) {
  assert(m_values != NULL);

  element.nodal_values(m_indices->array(), m_indices_e);
  for (int k = 0; k < element.n_chi(); k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      int i = 0, j = 0;
      element.local_to_global(k, i, j);
      x_nodal[k] = (*m_values)(i, j);
    }
  }
}

void DirichletData_Vector::enforce_homogeneous(const Element2 &element, Vector2d* x_nodal) {
  element.nodal_values(m_indices->array(), m_indices_e);
  for (int k = 0; k < element.n_chi(); k++) {
    if (m_indices_e[k] > 0.5) { // Dirichlet node
      x_nodal[k] = 0.0;
    }
  }
}

void DirichletData_Vector::fix_residual(Vector2d const *const *const x_global, Vector2d **r_global) {
  assert(m_values != NULL);

  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (auto p = grid.points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i] = m_weight * (x_global[j][i] - (*m_values)(i, j));
    }
  }
}

void DirichletData_Vector::fix_residual_homogeneous(Vector2d **r_global) {
  const IceGrid &grid = *m_indices->grid();

  // For each node that we own:
  for (auto p = grid.points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*m_indices)(i, j) > 0.5) {
      // Enforce explicit dirichlet data.
      r_global[j][i].u = 0.0;
      r_global[j][i].v = 0.0;
    }
  }
}

void DirichletData_Vector::fix_jacobian(Mat J) {
  const IceGrid &grid = *m_indices->grid();

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have taken steps to avoid touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.

  const double identity[4] = {m_weight, 0,
                              0, m_weight};
  ParallelSection loop(grid.com);
  try {
    for (auto p = grid.points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*m_indices)(i, j) > 0.5) {
        MatStencil row;
        row.j = j;
        row.i = i;
        PetscErrorCode ierr = MatSetValuesBlockedStencil(J, 1, &row, 1, &row, identity,
                                                         ADD_VALUES);
        PISM_CHK(ierr, "MatSetValuesBlockedStencil"); // this may throw
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

DirichletData_Vector::~DirichletData_Vector() {
  finish(m_values);
  m_values = NULL;
}

} // end of namespace fem
} // end of namespace pism
