/* Copyright (C) 2019, 2020, 2021, 2022 PISM Authors
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

#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "connected_components.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

/*!
 * Label connected components in a mask stored in an IceModelVec2S.
 *
 * This function allocates a copy on rank 0 and so should not be used if that is a
 * problem.
 */
void label_components(array::Scalar &mask, bool identify_icebergs, double mask_grounded) {
  auto mask_p0 = mask.allocate_proc0_copy();

  mask.put_on_proc0(*mask_p0);

  auto grid = mask.grid();

  ParallelSection rank0(grid->com);
  try {
    if (grid->rank() == 0) {
      petsc::VecArray array(*mask_p0);
      label_connected_components(array.get(),
                                 static_cast<int>(grid->My()),
                                 static_cast<int>(grid->Mx()),
                                 identify_icebergs,
                                 static_cast<int>(mask_grounded));
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  mask.get_from_proc0(*mask_p0);
}

} // end of namespace pism
