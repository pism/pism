/* Copyright (C) 2016 PISM Authors
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

#include "PISMAtmosphere.hh"

namespace pism {
namespace atmosphere {

AtmosphereModel::AtmosphereModel(IceGrid::ConstPtr g)
  : Component_TS(g),
    m_air_temp(m_sys, "effective_air_temp"),
    m_precipitation(m_sys, "effective_precipitation") {

  m_air_temp.set_string("pism_intent", "diagnostic");
  m_air_temp.set_string("long_name", "near-surface air temperature");
  m_air_temp.set_string("units", "K");

  m_precipitation.set_string("pism_intent", "diagnostic");
  m_precipitation.set_string("long_name", "precipitation");
  m_precipitation.set_string("units", "kg m-2 second-1");
  m_precipitation.set_string("glaciological_units", "kg m-2 year-1");
}

AtmosphereModel::~AtmosphereModel() {
  // empty
}

} // end of namespace atmosphere
} // end of namespace pism
