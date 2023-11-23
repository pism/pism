// Copyright (C) 2012-2014, 2016, 2023 PISM Authors
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

#pragma once

#include <memory>

#include "pism/hydrology/NullTransport.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace icebin {

class NullTransportHydrology : public pism::hydrology::NullTransport {
  friend class IBIceModel;

public:
  NullTransportHydrology(std::shared_ptr<const pism::Grid> grid);
  virtual ~NullTransportHydrology() = default;

  void update_impl(double icet, double icedt, const hydrology::Inputs& inputs);

protected:
  pism::array::Scalar basal_runoff_sum; // Cumulative effective thickness of water removed from till
};
} // end of namespace icebin
} // end of namespace pism
