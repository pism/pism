/* Copyright (C) 2015, 2020, 2023 PISM Authors
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

#ifndef _ICEGRID_REGIONAL_H_
#define _ICEGRID_REGIONAL_H_

namespace pism {

/**
 * Select a subset of the computational grid in an input file.
 *
 * @param ctx execution context
 *
 * @return A pointer to a new grid instance.
 */
std::shared_ptr<Grid> regional_grid_from_options(std::shared_ptr<Context> ctx);

} // end of namespace pism

#endif /* _ICEGRID_REGIONAL_H_ */
