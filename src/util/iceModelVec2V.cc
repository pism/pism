// Copyright (C) 2009--2017, 2020, 2021, 2022 Constantine Khroulev
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

#include "IceModelVec2V.hh"

#include "pism/util/array/Array_impl.hh"

#include "pism_utilities.hh"
#include "IceGrid.hh"

#include "pism/util/Context.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

IceModelVec2V::IceModelVec2V(IceGrid::ConstPtr grid, const std::string &name)
  : array::Array2D<Vector2>(grid, name, WITHOUT_GHOSTS, 2) {
  // This constructor uses the stencil width of 2 to make the DM compatible with ghosted
  // arrays with this wide stencil.

  auto sys = m_impl->grid->ctx()->unit_system();
  m_impl->metadata = {{sys, "u" + name}, {sys, "v" + name}};
  set_name("vel" + name);
}

IceModelVec2V::IceModelVec2V(IceGrid::ConstPtr grid, const std::string &name,
                             unsigned int stencil_width)
  : array::Array2D<Vector2>(grid, name,
                            stencil_width > 0 ? WITH_GHOSTS : WITHOUT_GHOSTS,
                            stencil_width) {

  auto sys = m_impl->grid->ctx()->unit_system();
  m_impl->metadata = {{sys, "u" + name}, {sys, "v" + name}};
  set_name("vel" + name);
}

std::shared_ptr<IceModelVec2V> IceModelVec2V::duplicate() const {

  auto result = std::make_shared<IceModelVec2V>(this->grid(),
                                                this->get_name());
  result->metadata(0) = this->metadata(0);
  result->metadata(1) = this->metadata(1);

  return result;
}

Velocity1::Velocity1(IceGrid::ConstPtr grid, const std::string &name)
  : IceModelVec2V(grid, name, 1) {
  // empty
}

Velocity1::Velocity1(IceGrid::ConstPtr grid, const std::string &name,
                     unsigned int stencil_width)
  : IceModelVec2V(grid, name, stencil_width) {
  // empty
}

Velocity2::Velocity2(IceGrid::ConstPtr grid, const std::string &name)
  : Velocity1(grid, name, 2) {
  // empty
}

} // end of namespace pism
