// Copyright (C) 2009--2017, 2020, 2021, 2022, 2023 Constantine Khroulev
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

#include "pism/util/array/Vector.hh"
#include "pism/util/array/Array_impl.hh"

#include "pism/util/Grid.hh"
#include "pism/util/Context.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {
namespace array {

Vector::Vector(std::shared_ptr<const Grid> grid, const std::string &name)
  : Array2D<pism::Vector2d>(grid, name, WITHOUT_GHOSTS, 2) {
  // This constructor uses the stencil width of 2 to make the DM compatible with ghosted
  // arrays with this wide stencil.

  auto sys = m_impl->grid->ctx()->unit_system();
  m_impl->metadata = {{sys, "u" + name}, {sys, "v" + name}};
  set_name("vel" + name);
}

Vector::Vector(std::shared_ptr<const Grid> grid, const std::string &name,
               unsigned int stencil_width)
  : Array2D<pism::Vector2d>(grid, name,
                           stencil_width > 0 ? WITH_GHOSTS : WITHOUT_GHOSTS,
                           stencil_width) {

  auto sys = m_impl->grid->ctx()->unit_system();
  m_impl->metadata = {{sys, "u" + name}, {sys, "v" + name}};
  set_name("vel" + name);
}

std::shared_ptr<Vector> Vector::duplicate() const {

  auto result = std::make_shared<Vector>(grid(), get_name());
  result->metadata(0) = this->metadata(0);
  result->metadata(1) = this->metadata(1);

  return result;
}

Vector1::Vector1(std::shared_ptr<const Grid> grid, const std::string &name)
  : Vector(grid, name, 1) {
  // empty
}

Vector1::Vector1(std::shared_ptr<const Grid> grid, const std::string &name,
                 unsigned int stencil_width)
  : Vector(grid, name, stencil_width) {
  // empty
}

Vector2::Vector2(std::shared_ptr<const Grid> grid, const std::string &name)
  : Vector1(grid, name, 2) {
  // empty
}

} // end of namespace array
} // end of namespace pism
