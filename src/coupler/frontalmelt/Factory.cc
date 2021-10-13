/* Copyright (C) 2015, 2017, 2018, 2019 PISM Authors
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

// frontal melt models:
#include "Constant.hh"
#include "DischargeGiven.hh"
#include "DischargeRouting.hh"
#include "Given.hh"

namespace pism {
namespace frontalmelt {
// FrontalMelt
Factory::Factory(IceGrid::ConstPtr g)
  : PCFactory<FrontalMelt>(g, "frontal_melt.models") {

  add_model<Constant>("constant");
  add_model<DischargeGiven>("discharge_given");
  add_model<DischargeRouting>("routing");
  add_model<Given>("given");
}

} // end of namespace frontalmelt
} // end of namespace pism
