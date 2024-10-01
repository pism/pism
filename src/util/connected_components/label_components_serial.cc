/* Copyright (C) 2019, 2020, 2021, 2022, 2023 PISM Authors
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

#include "label_components.hh"
#include "connected_components_impl.hh"

#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

namespace connected_components {

void label_serial(double *image, int nrows, int ncols, bool mark_isolated_components, int reachable,
                  int min_label) {

  using array = details::CArray<double>;

  array output(image, nrows, ncols);

  bool assign_final_labels = true;
  label(details::Mask<array>(output, reachable), mark_isolated_components, min_label,
        assign_final_labels, output);
}

/*!
 * Label connected components in a mask stored in an array::Scalar.
 *
 * This function allocates a copy on rank 0 and so should not be used if that is a
 * problem.
 */
void label_serial(array::Scalar &mask, bool mark_isolated_components, int reachable) {
  auto mask_p0 = mask.allocate_proc0_copy();

  mask.put_on_proc0(*mask_p0);

  auto grid = mask.grid();

  ParallelSection rank0(grid->com);
  try {
    if (grid->rank() == 0) {
      petsc::VecArray array(*mask_p0);
      int min_label = 1;
      connected_components::label_serial(array.get(), static_cast<int>(grid->My()),
                                         static_cast<int>(grid->Mx()), mark_isolated_components,
                                         reachable, min_label);
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  mask.get_from_proc0(*mask_p0);
}

} // namespace connected_components

} // end of namespace pism
