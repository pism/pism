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

#ifndef PISM_LABEL_COMPONENTS_H
#define PISM_LABEL_COMPONENTS_H

namespace pism {

namespace array {
class Scalar;
class Scalar1;
} // namespace array

namespace connected_components {

/*!
 * Label connected components in an `image`, modifying it "in place".
 *
 * The image is `nrows, ncols` in size. Positive values are treated as "foreground", zero
 * or negative as "background".
 *
 * If `mark_isolated_components` is set, then label isolated patches with `1` and patches
 * connected to values marked with `reachable` with zeros.
 *
 * If `mark_isolated_components` is *not* set, then assign a unique (consecutive) labels to each
 * patch, starting with `first_label`.
 */
void label_serial(double *image, int nrows, int ncols, bool mark_isolated_components, int reachable,
                  int min_label);

void label_serial(array::Scalar &mask, bool mark_isolated_components, int reachable);

void label(array::Scalar1 &mask);

void label_isolated(array::Scalar1 &mask, int reachable);

} // end of namespace connected_components

} // end of namespace pism

#endif /* PISM_LABEL_COMPONENTS_H */
