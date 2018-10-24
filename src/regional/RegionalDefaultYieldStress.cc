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

#include "RegionalDefaultYieldStress.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Vars.hh"

namespace pism {

RegionalDefaultYieldStress::RegionalDefaultYieldStress(IceGrid::ConstPtr g)
  : MohrCoulombYieldStress(g) {
}

RegionalDefaultYieldStress::~RegionalDefaultYieldStress() {
  // empty
}

void RegionalDefaultYieldStress::init_impl(const Geometry &geometry,
                                           const IceModelVec2S &till_water_thickness,
                                           const IceModelVec2S &overburden_pressure) {
  // turn off the second, redundant initialization message
  m_log->disable();
  MohrCoulombYieldStress::init_impl(geometry, till_water_thickness, overburden_pressure);
  m_log->enable();

  m_log->message(2,
                 "  using the regional version with strong till in no_model_mask area...\n");
}

void RegionalDefaultYieldStress::update_impl(const YieldStressInputs &inputs) {

  MohrCoulombYieldStress::update_impl(inputs);

  const IceModelVec2Int &nmm = *inputs.no_model_mask;

  double high_tauc = m_config->get_double("regional.no_model_yield_stress", "Pa");

  // now set tauc to a big value in no_model_strip
  IceModelVec::AccessList list{&nmm, &m_basal_yield_stress};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (nmm(i,j) > 0.5) {
      m_basal_yield_stress(i,j) = high_tauc;
    }
  }
}

} // end of namespace pism
