/* Copyright (C) 2015, 2016 PISM Authors
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

#include "RegionalDefaultYieldStress.hh"
#include "base/util/Logger.hh"
#include "base/util/PISMVars.hh"

namespace pism {

void RegionalDefaultYieldStress::init() {
  // turn off second, redundant init message
  m_log->disable();
  MohrCoulombYieldStress::init();
  m_log->enable();

  m_log->message(2,
             "  using the regional version with strong till in no_model_mask==1 area ...\n");
}

const IceModelVec2S& RegionalDefaultYieldStress::basal_material_yield_stress() {

  // do whatever you normally do
  const IceModelVec2S &result = MohrCoulombYieldStress::basal_material_yield_stress();

  // This is almost certainly redundant, but I don't want to count on
  // the fact that the base class puts results in m_basal_yield_stress.
  m_basal_yield_stress.copy_from(result);

  const IceModelVec2Int &nmm = *m_grid->variables().get_2d_mask("no_model_mask");

  // now set tauc to a big value in no_model_strip
  IceModelVec::AccessList list;
  list.add(nmm);
  list.add(m_basal_yield_stress);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (nmm(i,j) > 0.5) {
      m_basal_yield_stress(i,j) = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar
    }
  }
  return m_basal_yield_stress;
}

} // end of namespace pism
