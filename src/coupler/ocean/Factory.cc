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
#include "Anomaly.hh"
#include "Constant.hh"
#include "ConstantPIK.hh"
#include "GivenClimate.hh"
#include "Delta_T.hh"
#include "Delta_SMB.hh"
#include "Frac_MBP.hh"
#include "Frac_SMB.hh"
#include "Cache.hh"
#include "GivenTH.hh"
#include "Pico.hh"

namespace pism {
namespace ocean {
// Ocean
Factory::Factory(IceGrid::ConstPtr g)
  : PCFactory<OceanModel>(g) {
  m_option = "ocean";

  add_model<GivenTH>("th");
  add_model<PIK>("pik");
  add_model<Constant>("constant");
  add_model<Pico>("pico");
  add_model<Given>("given");
  set_default("constant");

  add_modifier<Anomaly>("anomaly");
  add_modifier<Cache>("cache");
  add_modifier<Delta_SMB>("delta_SMB");
  add_modifier<Frac_SMB>("frac_SMB");
  add_modifier<Delta_T>("delta_T");
  add_modifier<Frac_MBP>("frac_MBP");
}

Factory::~Factory() {
  // empty
}

} // end of namespace ocean
} // end of namespace pism
