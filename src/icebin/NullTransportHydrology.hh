// Copyright (C) 2012-2014, 2016 PISM Authors
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

#include <assert.h>

#include <base/hydrology/PISMHydrology.hh>
#include <base/stressbalance/PISMStressBalance.hh>
#include <base/util/PISMComponent.hh>
#include <base/util/iceModelVec.hh>
#include <base/util/iceModelVec2T.hh>

//! The PISM minimal model has till in a "can".  Water that overflows the can is not conserved.  There is no model for lateral transport.
/*!
This is the minimum functional derived class.  It updates till water thickness.

It has no transportable water and subglacial_water_thickness() returns zero.

This model can give no meaningful report on conservation errors.

Here is a talk which illustrates the "till-can" metaphor:
  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf
 */

namespace pism {
namespace icebin {

class NullTransportHydrology : public pism::hydrology::NullTransport {
  friend class IBIceModel;

public:
  NullTransportHydrology(pism::IceGrid::ConstPtr grid);
  virtual ~NullTransportHydrology() {
  }

  // solves an implicit step of a highly-simplified ODE
  void update_impl(double icet, double icedt);

protected:
  pism::IceModelVec2S basal_runoff_sum; // Cumulative effective thickness of water removed from till
};
} // end of namespace icebin
} // end of namespace pism
