/* Copyright (C) 2015 PISM Authors
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

#include "POFactory.hh"

// ocean models:
#include "POConstant.hh"
#include "POConstantPIK.hh"
#include "POGivenClimate.hh"
#include "PO_delta_SL.hh"
#include "PO_delta_T.hh"
#include "PO_delta_SMB.hh"
#include "PO_delta_MBP.hh"
#include "PO_frac_SMB.hh"
#include "POCache.hh"
#include "POGivenTH.hh" 

namespace pism {
namespace ocean {
// Ocean
Factory::Factory(IceGrid::ConstPtr g)
  : PCFactory<OceanModel,OceanModifier>(g) {
  m_option = "ocean";

  add_model<GivenTH>("th");
  add_model<PIK>("pik");
  add_model<Constant>("constant");
  add_model<Given>("given");
  set_default("constant");

  add_modifier<Cache>("cache");
  add_modifier<Delta_SMB>("delta_SMB");
  add_modifier<Frac_SMB>("frac_SMB");
  add_modifier<Delta_T>("delta_T");
  add_modifier<Delta_MBP>("delta_MBP");
  add_modifier<Delta_SL>("delta_SL");
}

Factory::~Factory() {
  // empty
}

} // end of namespace ocean
} // end of namespace pism
