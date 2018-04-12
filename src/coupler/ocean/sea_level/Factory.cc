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
#include "pism/coupler/SeaLevel.hh"
#include "Delta_SL.hh"
#include "Delta_SL_2D.hh"

namespace pism {
namespace ocean {
namespace sea_level {

Factory::Factory(IceGrid::ConstPtr grid)
  : PCFactory<SeaLevel>(grid) {
  m_option = "sea_level";

  add_model<SeaLevel>("constant");
  set_default("constant");

  add_modifier<Delta_SL>("delta_sl");
  add_modifier<Delta_SL_2D>("delta_sl_2d");
}

Factory::~Factory() {
  // empty
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
