/* Copyright (C) 2019, 2023 PISM Authors
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

#ifndef LABEL_COMPONENTS_H
#define LABEL_COMPONENTS_H

namespace pism {

namespace array {
class Scalar;
}

namespace petsc {
class Vec;
}

/*!
 * Label connected components in a mask stored in an array::Scalar.
 *
 * @param[in,out] mask_p0 temporary storage on rank 0.
 *
 * See the other `label_components()` details.
 */
void label_components(array::Scalar &mask,
                      petsc::Vec &mask_p0,
                      bool identify_icebergs, double mask_grounded);

/*!
 * Label connected components in a mask stored in an array::Scalar.
 *
 * This function allocates a copy on rank 0 and so should not be used if that is a
 * problem.
 *
 * @param[in,out] mask mask used to identify components (modified in place)

 * @param[in] identify_icebergs `true` to label blobs not connected to `mask_grounded` 1,
 *                              the rest with 0, `false` to assign unique labels to all blobs.
 * @param[in]  mask_grounded value in `mask` that is interpreted as "grounded"
 */
void label_components(array::Scalar &mask, bool identify_icebergs, double mask_grounded);

} // end of namespace pism

#endif /* LABEL_COMPONENTS_H */
