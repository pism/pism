/* Copyright (C) 2015, 2016, 2017, 2018 PISM Authors
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

#include "BedDef.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace bed {

Null::Null(IceGrid::ConstPtr g)
  : BedDef(g) {
  // empty
}

void Null::init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                     const IceModelVec2S &sea_level_elevation) {
  m_log->message(2,
             "* Initializing the dummy (no-op) bed deformation model...\n");

  BedDef::init_impl(opts, ice_thickness, sea_level_elevation);

  m_uplift.set(0.0);
}

MaxTimestep Null::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("bed_def none");
}

void Null::update_impl(const IceModelVec2S &ice_thickness,
                       const IceModelVec2S &sea_level_elevation,
                       double t, double dt) {
  (void) ice_thickness;
  (void) sea_level_elevation;
  (void) t;
  (void) dt;
  // This model does not update bed topography or bed uplift.
}

} // end of namespace bed
} // end of namespace pism
