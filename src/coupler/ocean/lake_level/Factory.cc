/* Copyright (C) 2015, 2017, 2018 PISM Authors
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

#include "Factory.hh"

// ocean models:
#include "pism/coupler/LakeLevel.hh"
#include "pism/coupler/ocean/lake_level/LakeCC.hh"

namespace pism {
namespace ocean {
namespace lake_level {

Factory::Factory(IceGrid::ConstPtr grid)
  : PCFactory<LakeLevel>(grid, "lake_level.models") {

  add_model<LakeLevel>("null");
  add_model<LakeCC>("lakecc");
}

Factory::~Factory() {
  // empty
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism