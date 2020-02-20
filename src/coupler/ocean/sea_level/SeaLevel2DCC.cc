/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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
#include <cassert>

#include "pism/util/pism_options.hh"
#include "pism/util/Time.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

#include "SeaLevel2DCC.hh"
#include "SeaLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace sea_level {

SeaLevel2DCC::SeaLevel2DCC(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in)
  : SeaLevel(g, in) {
}

SeaLevel2DCC::~SeaLevel2DCC() {
  // empty
}

void SeaLevel2DCC::init_impl(const Geometry &geometry) {
}

void SeaLevel2DCC::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

}

MaxTimestep SeaLevel2DCC::max_timestep_impl(double t) const {
}

bool SeaLevel2DCC::expandMargins_impl() const {
  return true;
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
