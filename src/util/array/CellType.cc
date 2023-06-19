/* Copyright (C) 2022, 2023 PISM Authors
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

#include "pism/util/array/CellType.hh"
#include "pism/util/interpolation.hh"

namespace pism {
namespace array {

CellType::CellType(std::shared_ptr<const Grid> grid, const std::string &name)
  : array::Scalar(grid, name) {
  set_interpolation_type(NEAREST);
}

CellType::CellType(std::shared_ptr<const Grid> grid, const std::string &name, int width)
  : array::Scalar(grid, name, width) {
  set_interpolation_type(NEAREST);
}

CellType1::CellType1(std::shared_ptr<const Grid> grid, const std::string &name)
  : CellType(grid, name, 1) {
  // empty
}

CellType1::CellType1(std::shared_ptr<const Grid> grid, const std::string &name, int width)
  : CellType(grid, name, width) {
  // empty
}

CellType2::CellType2(std::shared_ptr<const Grid> grid, const std::string &name)
  : CellType1(grid, name, 2) {
  // empty
}

} // end of namespace array
} // end of namespace pism
