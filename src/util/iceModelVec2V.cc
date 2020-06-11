// Copyright (C) 2009--2017, 2020 Constantine Khroulev
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

#include <memory>               // std::dynamic_pointer_cast

#include "IceModelVec2V.hh"
#include "IceModelVec_impl.hh"

#include "pism_utilities.hh"
#include "IceGrid.hh"

#include "error_handling.hh"    // RuntimeError

#include "iceModelVec_helpers.hh" // add_2d, copy_2d

#include "pism/util/Context.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

IceModelVec2V::IceModelVec2V(IceGrid::ConstPtr grid, const std::string &short_name,
                             IceModelVecKind ghostedp, unsigned int stencil_width)
  : IceModelVec2Struct<Vector2>(grid, short_name, ghostedp, stencil_width) {

  auto sys = m_impl->grid->ctx()->unit_system();

  m_impl->metadata = {SpatialVariableMetadata(sys, "u" + short_name),
                      SpatialVariableMetadata(sys, "v" + short_name)};

  m_impl->name = "vel" + short_name;
}

IceModelVec2V::~IceModelVec2V() {
  // empty
}

IceModelVec2V::Ptr IceModelVec2V::ToVector(IceModelVec::Ptr input) {

 IceModelVec2V::Ptr result = std::dynamic_pointer_cast<IceModelVec2V,IceModelVec>(input);

  if (not (bool)result) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
  }

  return result;
}

void IceModelVec2V::add(double alpha, const IceModelVec &x) {
  return add_2d<IceModelVec2V>(this, alpha, &x, this);
}

void IceModelVec2V::add(double alpha, const IceModelVec &x, IceModelVec &result) const {
  return add_2d<IceModelVec2V>(this, alpha, &x, &result);
}

void IceModelVec2V::copy_from(const IceModelVec &source) {
  return copy_2d<IceModelVec2V>(&source, this);
}

} // end of namespace pism
