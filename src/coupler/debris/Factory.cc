/* Copyright (C) 2026 PISM Authors
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

#include "pism/coupler/debris/Factory.hh"

// debris models:
#include "pism/coupler/debris/GivenDebris.hh"

namespace pism {
namespace debris {

Factory::Factory(std::shared_ptr<const Grid> g)
  : PCFactory<DebrisModel>(g, "debris.models") {

  add_model<Given>("given");
}

} // end of namespace debris
} // end of namespace pism
