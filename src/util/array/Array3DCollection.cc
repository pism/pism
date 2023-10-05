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
#include <vector>

#include "pism/util/array/Array3DCollection.hh"

#include "pism/util/Grid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/array/Array_impl.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Context.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/LocalInterpCtx.hh"
#include "pism/util/io/io_helpers.hh"

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

void Array3DCollection::regrid_impl(const File &file, io::Default default_value) {

  auto log = grid()->ctx()->log();

  auto variable = metadata(0);

  auto V = file.find_variable(variable.get_name(), variable["standard_name"]);

  petsc::TemporaryGlobalVec tmp(dm());

  if (not V.exists) {
    set_default_value_or_stop(file.filename(), variable, default_value, *log, tmp);
  } else {
    bool allow_extrapolation = grid()->ctx()->config()->get_flag("grid.allow_extrapolation");

    grid::InputGridInfo input_grid(file, V.name, variable.unit_system(), grid()->registration());

    // The lengths of third dimensions of a "collection" in the input file and internally
    // *have to match*.
    {
      auto dim_name          = variable.z().get_name();
      auto n_levels_internal = get_levels().size();
      auto n_levels_input    = file.dimension_length(dim_name);

      if (n_levels_input != n_levels_internal) {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "the length of the '%s' dimension in '%s' has to match the length of '%s' used internally"
            " (got %d and %d)",
            dim_name.c_str(), file.filename().c_str(), dim_name.c_str(), (int)n_levels_input,
            (int)n_levels_internal);
      }
    }

    input_grid.report(*log, 4, variable.unit_system());

    // Note: Array3DCollection is a collection of 2D fields, so we don't need to use the
    // internal vertical (Z) grid.
    io::check_input_grid(input_grid, *grid(), {0.0});

    // initialize interpolation in X and Y directions:
    LocalInterpCtx lic(input_grid, *grid(), m_impl->interpolation_type);

    // fake the Z axis: the values of the corresponding coordinate variable in the file
    // may not be unique or increasing.
    std::vector<double> Z;
    {
      auto name = variable.z().get_name();
      auto N = file.dimension_length(name);
      Z.resize(N);
      for (int k = 0; k < N; ++k) {
        Z[k] = k;
      }
    }
    lic.start[Z_AXIS] = 0;
    lic.count[Z_AXIS] = Z.size();
    lic.z = std::make_shared<Interpolation>(NEAREST, Z, Z);

    // Note: this call will read the last time record (the index is set in `lic` based on
    // info in `input_grid`).
    petsc::VecArray tmp_array(tmp);
    io::regrid_spatial_variable(variable, input_grid, *grid(), lic, file, tmp_array.get());
  }

  if (m_impl->ghosted) {
    global_to_local(*dm(), tmp, vec());
  } else {
    PetscErrorCode ierr = VecCopy(tmp, vec());
    PISM_CHK(ierr, "VecCopy");
  }
}

} // end of namespace array
} // end of namespace pism
