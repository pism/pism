// Copyright (C) 2008--2018, 2020, 2021, 2022, 2023 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <memory>

#include "pism/util/array/Array3DCollection.hh"

#include "pism/util/Grid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/array/Array_impl.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace array {

// this file contains method for derived class array::Array3DCollection

Array3DCollection::Array3DCollection(std::shared_ptr<const Grid> grid, const std::string &name,
                                     Kind ghostedp, const std::vector<double> &levels,
                                     unsigned int stencil_width)
    : Array(grid, name, ghostedp, 1, stencil_width, levels) {
  set_begin_access_use_dof(true);
}

double *Array3DCollection::column(int i, int j) {
#if (Pism_DEBUG == 1)
  check_array_indices(i, j, 0);
#endif
  return ((double ***)m_array)[j][i];
}

const double *Array3DCollection::column(int i, int j) const {
#if (Pism_DEBUG == 1)
  check_array_indices(i, j, 0);
#endif
  return ((double ***)m_array)[j][i];
}

void Array3DCollection::copy_from(const Array3DCollection &input) {
  array::AccessScope list{ this, &input };

  assert(get_levels().size() == input.get_levels().size());
  assert(ndof() == input.ndof());

  auto N = get_levels().size();

  ParallelSection loop(m_impl->grid->com);
  try {
    for (auto p = m_impl->grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

#if PETSC_VERSION_LT(3, 12, 0)
      PetscMemmove(this->column(i, j), const_cast<double *>(input.column(i, j)),
                   N * sizeof(double));
#else
      PetscArraymove(this->column(i, j), input.column(i, j), N);
#endif
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  update_ghosts();

  inc_state_counter();
}

std::shared_ptr<Array3DCollection> Array3DCollection::duplicate() const {

  auto result = std::make_shared<Array3DCollection>(this->grid(), this->get_name(), WITHOUT_GHOSTS,
                                                    this->get_levels());

  result->metadata() = this->metadata();

  return result;
}

} // end of namespace array
} // end of namespace pism
