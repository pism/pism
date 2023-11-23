/* Copyright (C) 2016, 2019 PISM Authors
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

#ifndef REMOVE_NARROW_TONGUES_H
#define REMOVE_NARROW_TONGUES_H

namespace pism {

namespace array { class Scalar; }
class Geometry;

void remove_narrow_tongues(const Geometry &geometry, array::Scalar &ice_thickness);

} // end of namespace pism


#endif /* REMOVE_NARROW_TONGUES_H */
