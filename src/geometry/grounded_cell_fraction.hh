/* Copyright (C) 2016, 2017, 2018, 2022 PISM Authors
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

#ifndef _GROUNDED_CELL_FRACTION_H_
#define _GROUNDED_CELL_FRACTION_H_

namespace pism {

namespace array { class Scalar; }

double grounded_area_fraction(double a, double b, double c);

/*!
 * Compute grounded cell fractions by splitting control volumes into triangles and
 * treating the flotation criterion as a linear function on each triangle.
 */
void compute_grounded_cell_fraction(double ice_density,
                                    double ocean_density,
                                    const array::Scalar &sea_level,
                                    const array::Scalar &ice_thickness,
                                    const array::Scalar &bed_topography,
                                    array::Scalar &result);

} // end of namespace pism


#endif /* _GROUNDED_CELL_FRACTION_H_ */
